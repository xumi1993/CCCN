#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sacio.h"

void filter4_(double *f1,double *f2,double *f3,double *f4,int *npow,
              double *dt,int *n, float seis_in[],float seis_out[],
              float seis_outamp[],
              float seis_outph[],int *ns,double *dom);

int main(int argc, char **argv){
    static int n, ns,npow;
    static double f1, f2, f3, f4, dt,dom;
    static float seis_in[10000000],seis_out[10000000];
    static float seis_outamp[10000000],seis_outph[10000000];
    SACHEAD hd;
    char  name[160];
    float *trace, *out_real, *out_imag, *out_trace;
    int i;

    strcpy(name,argv[1]);
    trace = read_sac(name, &hd);
    n  = hd.npts;
    dt = hd.delta;
    npow = 1;
    for (i=0;i<n;i++){
        seis_in[i] = trace[i];
    }
    f1 = 0.01;f2=0.02;f3=0.067;f4=0.08;
    filter4_(&f1,&f2,&f3,&f4,&npow,&dt,&n,seis_in,seis_out,seis_outamp,seis_outph,&ns,&dom);
    hd = new_sac_head(dom,ns,0);
    out_real = (float*) calloc(ns,sizeof(float));
    out_imag = (float*) calloc(ns,sizeof(float));
    for (i=0;i<ns;i++){
        out_real[i] = seis_outamp[i]*cos(seis_outph[i]);
        out_imag[i] = seis_outamp[i]*sin(seis_outph[i]);
    }
    write_sac("out.rl",hd,out_real);
    write_sac("out.im",hd,out_imag);
    free(out_real);
    free(out_imag);
    
    out_trace = (float*) calloc(n,sizeof(float));
    hd = new_sac_head(dt,n,0);
    for (i=0;i<n;i++){
        out_trace[i] = seis_out[i];
    }
    write_sac("out.tr",hd,out_trace);
    free(out_trace);
}
