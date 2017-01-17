#define MAIN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mysac.h"

/* Function prorotypes */


void filter4_(double *f1,double *f2,double *f3,double *f4,int *npow,
              double *dt,int *n, float seis_in[],float seis_out[],
              float seis_outamp[],
              float seis_outph[],int *ns,double *dom);

void swapn(unsigned char *b, int N, int n);




/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	SAC_HD *read_sac (char *fname, float *sig, SAC_HD *SHD, long nmax)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
/*..........................................................................*/
	fsac = fopen(fname, "rb");
	if ( !fsac )
	{
	 fclose (fsac);
	 return NULL;
	}

	if ( !SHD ) SHD = &SAC_HEADER;

	 fread(SHD,sizeof(SAC_HD),1,fsac);

	 if ( SHD->npts > nmax )
	 {
	  fprintf(stderr,
	   "ATTENTION !!! dans le fichier %s npts est limite a %d",fname,nmax);

	  SHD->npts = nmax;
	 }

	 fread(sig,sizeof(float),(int)(SHD->npts),fsac);

	fclose (fsac);

   /*-------------  calcule de t0  ----------------*/
   {
	int eh, em ,i;
	float fes;
	char koo[9];

	for ( i = 0; i < 8; i++ ) koo[i] = SHD->ko[i];
	koo[8] = NULL;

	SHD->o = SHD->b + SHD->nzhour*3600. + SHD->nzmin*60 +
	 SHD->nzsec + SHD->nzmsec*.001;

	sscanf(koo,"%d%*[^0123456789]%d%*[^.0123456789]%g",&eh,&em,&fes);

	SHD->o  -= (eh*3600. + em*60. + fes);
   /*-------------------------------------------*/}

	return SHD;
}

/*c/////////////////////////////////////////////////////////////////////////*/
/*--------------------------------------------------------------------------*/
	void write_sac (char *fname, float *sig, SAC_HD *SHD)
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
{
 FILE *fsac;
 int i;
/*..........................................................................*/
	fsac = fopen(fname, "wb");

	if ( !SHD ) SHD = &SAC_HEADER;


        SHD->iftype = (long)ITIME;
        SHD->leven = (long)TRUE;

        SHD->lovrok = (long)TRUE;
        SHD->internal4 = 6L;



  /*+++++++++++++++++++++++++++++++++++++++++*/
     SHD->depmin = sig[0];
     SHD->depmax = sig[0];
 
   for ( i = 0; i < SHD->npts ; i++ )
   {
    if ( SHD->depmin > sig[i] ) SHD->depmin = sig[i];
    if ( SHD->depmax < sig[i] ) SHD->depmax = sig[i];
   }

	 fwrite(SHD,sizeof(SAC_HD),1,fsac);

	 fwrite(sig,sizeof(float),(int)(SHD->npts),fsac);


	 fclose (fsac);
}



/*c/////////////////////////////////////////////////////////////////////////*/
 float sig[10000000];
 SAC_HD shd1;


/*c/////////////////////////////////////////////////////////////////////////*/

int main (int argc, char *argv[])
{
static int n, ns,npow;
static double f1, f2, f3, f4, dt,dom;
static float seis_in[10000000],seis_out[10000000];
static float seis_outamp[10000000],seis_outph[10000000];
double t1,t2,t3,t4;
char  name[160], name1[160];
char  nameamp[160],nameph[160];
FILE  *in, *ff;
int   i, j, nn;


  if( argc != 2) {
      printf("Usage: whiten_phamp  parameter_file\n");
      exit(-1);
  }

// open and read parameter file param.dat
  if((in = fopen(argv[1],"r")) == NULL) {
      printf("Can not find file %s.\n",argv[1]);
      exit(1);
  }

  while((nn = fscanf(in,"%lf %lf %lf %lf %lf %d %s",&t1,&t2,&t3,&t4,
            &dt,&npow,name)) != EOF) { // start main loop
      if(nn == 0 || nn != 7) break;
      printf("Corners periods. Low: %f - %f, High: %f - %f\n",t1, t2, t3, t4);
      printf("Step: %f, Cosine power: %d\n",dt, npow);

// remove quotes from name
      j = 0;
      for(i = 0; i < strlen(name); i++) {
          if(name[i] == '\'' || name[i] == '\"') continue;
          name[j] = name[i]; j++;
      }
      name[j] = '\0';

// do running average before whitening

  ff = fopen("sac_one_cor","w");
//  fprintf(ff,"/home/nshapiro/PROGS/SAC/bin/sac2000 << END\n");
   fprintf(ff,"sac << END\n");
  fprintf(ff,"r %s\n", name);
  fprintf(ff,"abs\n");
  fprintf(ff,"smooth mean h 40\n");
  fprintf(ff,"w a.avg \n");
  fprintf(ff,"r %s\n",name);
  fprintf(ff,"divf a.avg\n");
  fprintf(ff,"w smooth.sac\n");
  fprintf(ff,"quit\n");
  fprintf(ff,"END\n");
  fclose(ff);
  system("sh sac_one_cor");

  // end of running average


  if ( !read_sac("smooth.sac", sig, &shd1, 10000000 ) )
    {
      fprintf(stderr,"file %s did not found\n", name1 );
      return 0;
    }

    n  = shd1.npts;
    dt = shd1.delta;
   
     for( i =0; i< n; i++)
     {  
     seis_in[i] = sig[i];  
     //     printf(" seis_in1  %d %f\n", i,sig[i]);
     }

      printf(" Dt1= %f, Nsamples1= %d\n",dt, n);


  f1 = 1.0/t1; f2 = 1.0/t2; f3 = 1.0/t3; f4 = 1.0/t4;
  filter4_(&f1,&f2,&f3,&f4,&npow,&dt,&n,seis_in,seis_out,seis_outamp,seis_outph,&ns,&dom);
     

        write_sac(name,seis_out, &shd1);

        strcpy(nameamp,name);
        strcpy(nameph,name);
        strcat(nameamp,".am");
        strcat(nameph, ".ph");
	shd1.npts = ns/2 + 1;
	shd1.delta = dom;
	shd1.b = 0;
        shd1.iftype = IXY;
        write_sac(nameamp,seis_outamp, &shd1 );
        write_sac(nameph, seis_outph,  &shd1 );



  }
  
   return 0;
  }
