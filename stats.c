#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <getopt.h>

int main(int argc, char **argv)
{
  double from=0;
  double to=1000;
  double step=0.1;
  long steps;
  int flags, opt;
  int minlen=500;
  
  while ((opt = getopt(argc, argv, "f:t:s:l:")) != -1) {
    switch (opt) {
    case 'f':
      from=atof(optarg);
      break;
    case 't':
      to=atof(optarg);
      break;
    case 's':
      step = atof(optarg);
      break;
    case 'l':  
      minlen=atoi(optarg);
      break;
      
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-f from time] [-t to time] [-s time step] [-l min length]\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  steps=(long)((to-from)/step)+1;
  long *allvisits=(long *)calloc(steps,sizeof(long));
  double *allxin=(double*)calloc(steps,sizeof(double));
  double *allyn=(double*)calloc(steps,sizeof(double));
  double *allxinsquare=(double*)calloc(steps,sizeof(double));
  double *allynsquare=(double*)calloc(steps,sizeof(double));
  double *allxinyn=(double*)calloc(steps,sizeof(double));

  long *visits;
  double *xin;
  double *yn;
  double *xinsquare;
  double *ynsquare;
  double *xinyn;
  long len;
  int err;

  double t,xi,y;
  long index,i;

  err=scanf("{");
  while(err!=EOF){
    err=scanf("{");


    visits=(long *)calloc(steps,sizeof(long));
    xin=(double*)calloc(steps,sizeof(double));
    yn=(double*)calloc(steps,sizeof(double));
    xinsquare=(double*)calloc(steps,sizeof(double));
    ynsquare=(double*)calloc(steps,sizeof(double));
    xinyn=(double*)calloc(steps,sizeof(double));
    len=0;
    while((err=scanf("{%lf,%lf,%lf},",&t,&xi,&y))==3) {

      index=(long)(floor((t-from)/step));      
      //      fprintf(stderr,"%ld %.10g %.10g %.10g\n",index,t,xi,y);      
      if ((index>=0) && (index<steps)) {
	visits[index]++;
	xin[index]+=xi;
	yn[index]+=y;
	xinsquare[index]+=xi*xi;
	ynsquare[index]+=y*y;
	xinyn[index]+=xi*y;
      }
      len++;
    }
    if (len>=minlen) {
      for (i=0;i<steps;i++) {
	if (visits[i]>0) {
	  allvisits[i]++;
	  allxin[i]+=xin[i]/visits[i];
	  allyn[i]+=yn[i]/visits[i];
	  allxinsquare[i]+=xinsquare[i]/visits[i]; //xinsquare[i]/visits[i];
	  allynsquare[i]+=ynsquare[i]/visits[i];//ynsquare[i]/visits[i];
	  allxinyn[i]+=xinyn[i]/visits[i];
	}
      }
    }
    free(visits);
    free(xin);
    free(yn);
    free(xinsquare);
    free(ynsquare);
    free(xinyn);
    
    err=scanf("}, ");
  }

  for (i=0;i<steps;i++) {
    if (allvisits[i]>0) {
      printf("%.10g %ld %.10g %.10g %.10g %.10g %.10g\n",from+step*i,allvisits[i],allxin[i]/allvisits[i],allxinsquare[i]/allvisits[i],allyn[i]/allvisits[i],allynsquare[i]/allvisits[i],allxinyn[i]/allvisits[i]);
    }
  }

  free(allvisits);
  free(allxin);
  free(allyn);
  free(allxinsquare);
  free(allynsquare);
  free(allxinyn);
  
  return 0;
}  
