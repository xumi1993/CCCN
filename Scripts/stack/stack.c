/*************************************************/
/*FileName: stack.c                              */
/*Author  : xfeng                                */
/*Mail    : geophydogvon@gmail.com               */
/*Inst    : NJU                                  */
/*Time    :Tue 14 Feb 2017 03:05:02 PM CST       */
/*This is c programming language!                */
/*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sacio.h"


void str_no_void(char *ps){
    char *pt = ps;
    while ( *ps != '\0' ){
        if( *ps != ' ' && *ps != '\n' ){
            *pt++ = *ps;
        }
        ++ps;
    }
    *pt = '\0';
}


int main(int argc, char *argv[]){
    int i = 0, j = 0, size = 256;
    float *data, sum[500000] = {0.},*outdata, starttime, endtime;
    char* ss = (char*)malloc(size);
    FILE *fp;
    SACHEAD hd;

    starttime = clock();

    if( argc != 3 ){
        fprintf(stderr,"Usage: stack <sacfile_lsit> <out_stacked_sac>\n");
        fprintf(stderr,"       return sacfile after stacking inputing files\n");
        exit(1);
    }

    if( (fp = fopen(argv[1],"r")) == NULL ){
        fprintf(stderr,"Error reading file: %s\n",argv[1]);
        exit(1);
    }
    while( fgets(ss,size,fp) ){
        str_no_void(ss);
        //printf("Reading file ->>>>>>>>>>> %s\n",ss);
        data = read_sac(ss,&hd);
        for ( i = 0; i < hd.npts; i ++ ){
            sum[i] += data[i];
        }
        free(data);
    }

    outdata = (float*) malloc(sizeof(float)*hd.npts);
    for ( i = 0; i < hd.npts; i ++ ){
        outdata[i] = sum[i];
    }
    write_sac(argv[2],hd,outdata);
    fclose(fp);
    free(outdata);
    endtime = clock();
    printf("Stacked out file %s costs %f seconds!\n",argv[2],(endtime-starttime)/CLOCKS_PER_SEC);
    return 0;
}
