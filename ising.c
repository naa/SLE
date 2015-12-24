#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <gmp.h>

gmp_randstate_t state;
mpf_t rop;

int **allocate_2d_array(int size)
{
  int **table;
  table = (int**)malloc(size * sizeof(int*));
  for (int i = 0; i < size; i++) {
    table[i] = (int*)malloc(size * sizeof(int));
  }
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      table[i][j] = 0;
    }
  }
  return table;
}

void free_2d_array(int size, int **table){
  for (int i = 0; i < size; i++) {
    free(table[i]);
  }
  free(table);
}

int mod(int x, int m) {
    return (x%m + m)%m;
}

int rand_int(int max) 
{  
  return (int)floor(((double) max * gmp_urandomm_ui(state,(long)RAND_MAX+1)) / (RAND_MAX+1.0));
}

double nrand()
{
  mpf_urandomb (rop,  state, 256);
  return mpf_get_d(rop);
  //  return (double)(gmp_urandomm_ui(state,(long)RAND_MAX+1)) / (double)RAND_MAX ;
}

void init_table_pm_boundary(int size, int ** table)
{
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) { 
      table[i][j] = 2 * (nrand() > 0.5) - 1;
      if (( i == 0 ) || ( ( j == size-1 ) && ( i < size - 1)))
	table[i][j]=-1;
      if (( i == size-1 ) || (( j == 0 ) && ( i > 0 )))
	table[i][j]=1;
    }
  }
}

void init_table_no_boundary(int size, int ** table)
{
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) { 
      table[i][j] = 2 * (nrand() > 0.5) - 1;
      //table[i][j] = 1;//-1+2*((i+j)%2);
    }
  }
}

int modify_cell (int size,int **table, int i, int j, double beta)
{
  if (table[mod(i,size)][mod(j,size)] ==0) return 0;
  double p1 = exp( beta * (table[mod(i - 1,size)][mod(j,size)] + 
			   table[mod(i + 1,size)][mod(j,size)] + 
			   table[mod(i,size)][mod(j - 1,size)] + 
			   table[mod(i,size)][mod(j + 1,size)]));
  double p2 = exp(-beta * (table[mod(i - 1,size)][mod(j,size)] + 
			   table[mod(i + 1,size)][mod(j,size)] + 
			   table[mod(i,size)][mod(j - 1,size)] + 
			   table[mod(i,size)][mod(j + 1,size)]));
  int nval = 2 * (((p1 + p2) * nrand()) < p1) - 1;
  int res = (nval != table[mod(i,size)][mod(j,size)]);
  table[mod(i,size)][mod(j,size)]=nval;
  return res;
}

//int modify_cell (int size,int **table, int i, int j, double beta)
//{
//  if (table[mod(i,size)][mod(j,size)]==0) return 0;
//  double s1,s2;
//  s1=(table[mod(i - 1,size)][mod(j,size)] + 
//      table[mod(i + 1,size)][mod(j,size)] + 
//      table[mod(i,size)][mod(j - 1,size)] + 
//      table[mod(i,size)][mod(j + 1,size)])*
//    table[mod(i,size)][mod(j,size)];
//  s2=-s1;
//  if (((s2-s1)<0) || (nrand()<exp(-beta*(s2-s1))))
//    table[mod(i,size)][mod(j,size)]=-table[mod(i,size)][mod(j,size)];
//  return 0;
//}
//
void modify_cluster (int size,int **table, int i, int j, double beta)
{
  if (table[mod(i,size)][mod(j,size)]==0) return;   
  int *qx=(int*)malloc(size * size * sizeof(int));
  int *qy=(int*)malloc(size * size * sizeof(int));  
  long head=0;					
  long tail=0;
  int x,y;
  long blength=0;
  int k;
  int cluster=table[mod(i,size)][mod(j,size)];
  table[mod(i,size)][mod(j,size)]=3;
  qx[tail]=i;
  qy[tail]=j;
  tail++;
  while (tail>head) {
    x=qx[head];
    y=qy[head];
    head++;
    if ((table[mod(x-1,size)][mod(y,size)]==cluster) && (nrand()<1.0-exp(-2*beta))){
      qx[tail]=mod(x-1,size);
      qy[tail]=y;
      table[mod(x-1,size)][mod(y,size)]=3;
      tail++;
    } else if (table[mod(x-1,size)][mod(y,size)]==-cluster) {
      blength++;
    }
    if ((table[mod(x+1,size)][mod(y,size)]==cluster) && (nrand()<1.0-exp(-2*beta))){
      qx[tail]=mod(x+1,size);
      qy[tail]=y;
      table[mod(x+1,size)][mod(y,size)]=3;
      tail++;
    } else if (table[mod(x+1,size)][mod(y,size)]==-cluster) {
      blength++;
    }
    if ((table[mod(x,size)][mod(y-1,size)]==cluster) && (nrand()<1.0-exp(-2*beta))){
      qx[tail]=mod(x,size);
      qy[tail]=mod(y-1,size);
      table[mod(x,size)][mod(y-1,size)]=3;      
      tail++;
    } else if (table[mod(x,size)][mod(y-1,size)]==-cluster) {
      blength++;
    }
    if ((table[mod(x,size)][mod(y+1,size)]==cluster) && (nrand()<1.0-exp(-2*beta))){
      qx[tail]=mod(x,size);
      qy[tail]=mod(y+1,size);
      table[mod(x,size)][mod(y+1,size)]=3;            
      tail++;
    } else if (table[mod(x,size)][mod(y+1,size)]==-cluster) {
      blength++;
    }
  }
  double p1 = exp( beta * blength);
  double p2 = exp(-beta * blength);
  //int nval = 2 * (((p1 + p2) * nrand()) < p1) - 1;
  int nval = -cluster;
  //  printf("\n%lf %d %d %ld %ld %ld\n", 1/beta, cluster, nval, blength, head, tail);  
  for (k=0; k<head; k++) {
    table[qx[k]][qy[k]]=nval;
  }
  free(qx);
  free(qy);
}

void evolve_table_cyclic_boundary ( int size, int** table, double beta, long long niter ) 
{
  int i,j;

  for (long long k = 0; k < niter; k++){
    for (i=0;i<size;i++)
      for (j=0;j<size;j++)
	//modify_cluster (size,table,rand_int(size-1), rand_int(size-1),beta);
	modify_cell (size,table,rand_int(size), rand_int(size),beta);
	//modify_cell (size,table,i,j,beta);
  }
}

void evolve_table_pm_boundary ( int size, int** table, double beta, long long niter ) 
{
  int i,j;
  init_table_pm_boundary(size,table);

  for (long long k = 0; k < niter; k++){
    i=1 + rand_int(size-2);
    j=1 + rand_int(size-2);
    modify_cell (size,table,i, j,beta);
  }
}

int **  calc_interface(int size, int **table, int ** hor_interface, int ** vert_interface)
{
  int **interface = allocate_2d_array(size);

  int x1, y1, x2, y2, dx, dy; // dx = x2 - x1; dy = y2 - y1;
  int c1, c2;

  for (x1 = 0; (x1 < size - 1) && (table[x1][0] == table[x1 + 1][0]); x1++) {}
  x2 = x1 + 1; 
  y1 = 0; 
  y2 = 0; 

  for (int k = 0; ((y2 < size - 1) || (y1 < size - 1) ) && (k < size * size); k++)    {
    interface[x1][y1] = 1;
    interface[x2][y2] = 2;
    if (x1==x2) vert_interface[x1][y1]=1;
    if (y1==y2) hor_interface[x1][y1]=1;

    dx = x2 - x1; dy = y2 - y1;
    c1 = table[x1 - dy][y1 + dx];
    c2 = table[x2 - dy][y2 + dx];
    if ((3 * c1 + c2) == (3 * (+1) + (+1))) {x2 = x1 - dy; y2 = y1 + dx;}
    if ((3 * c1 + c2) == (3 * (-1) + (-1))) {x1 = x2 - dy; y1 = y2 + dx;}
    if ((3 * c1 + c2) == (3 * (-1) + (+1))) {x1 = x1 - dy; y1 = y1 + dx; x2 = x2 - dy; y2 = y2 + dx;}
    if ((3 * c1 + c2) == (3 * (+1) + (-1))) {x1 = x2 - dy; y1 = y2 + dx;}
  }
  interface[x1][y1] = 1;
  interface[x2][y2] = 2;
  if (x1==x2) vert_interface[x1][y1]=1;
  if (y1==y2) hor_interface[x1][y1]=1;

  return interface;
}

void add_crossing_number(int size, int ** hor_interface, int ** res_table) 
{
  int num=0;
  for (int j=0; j<size; j++) {
    num=0;
    for (int i=0; i<size; i++) {
      res_table[i][j] += num%2;
      if (hor_interface[i][j]!=0) 
	num++;
    }
  }
}

long calc_magnetization(int size, int ** table)
{
  long res=0;
  for (int i=0; i<size; i++)
    for (int j=0; j<size; j++)
      res+=table[i][j];
  return res;
}

void print_table (FILE* F, int size, int **table)
{
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (table[i][j] > 0) 
	fprintf(F, "+"); 
      else
	fprintf(F, "-");
    }
    fprintf(F, "\n");
  }
}

void print_table_xpm_simple (FILE* F, int size, int **table)
{
  fprintf(F, "/* XPM */\nstatic char *pict[] = {\n");
  fprintf(F, "\"%d %d %d %d\",\n", size, size, 2, 1);
  fprintf(F, "\"+ c #ffffff\",\n");
  fprintf(F, "\"- c #000000\",\n");

  for (int i = 0; i < size; i++) {
    fprintf(F, "\"");
    for (int j = 0; j < size; j++)
      if (table[i][j] > 0) 
	fprintf(F, "+"); 
      else 
	fprintf(F, "-");
    fprintf(F, "\"");
    if (i < size - 1) {
      fprintf(F, ",");
    }
    fprintf(F, "\n");
  }

  fprintf(F, "};");
}


void print_table_xpm (FILE* F, int size, int **table, int **interface)
{
  fprintf(F, "/* XPM */\nstatic char *pict[] = {\n");
  fprintf(F, "\"%d %d %d %d\",\n", size, size, 4, 1);
  fprintf(F, "\"+ c #ffffff\",\n");
  fprintf(F, "\"- c #000000\",\n");
  //fprintf(F, "\"2 c #ffcfcf\",\n");
  //fprintf(F, "\"1 c #300000\",\n");
  fprintf(F, "\"2 c #ff0000\",\n");
  fprintf(F, "\"1 c #ff0000\",\n");

  for (int i = 0; i < size; i++) {
    fprintf(F, "\"");
    for (int j = 0; j < size; j++)
      if (interface[i][j]==0) {
	if (table[i][j] > 0) 
	  fprintf(F, "+"); 
	else 
	  fprintf(F, "-");
      } else {
	fprintf(F, "%d", interface[i][j]);
      }
    fprintf(F, "\"");
    if (i < size - 1) {
      fprintf(F, ",");
    }
    fprintf(F, "\n");
  }

  fprintf(F, "};");
}

void print_table_as_list (FILE* F, int size,long measurements, int **table)
{
  fprintf(F, "{");
  for (int i = 0; i < size; i++) {
    fprintf(F,"{");
    for (int j = 0; j < size; j++) {
      fprintf(F, "%lf",((double)table[i][j])/measurements); 
      if (j<size-1)
	fprintf(F, ",");
    }
    fprintf(F, "}");
    if (i<size-1)
      fprintf(F,",");
    fprintf(F, "\n");
  }
  fprintf(F, "}");
}

void compute_probability(int size, double beta, long measurements,int **res) 
{
  long long niter=3*size*size;
  int **table=allocate_2d_array(size);
  //  int ** res = allocate_2d_array(size);
  int ** hint=allocate_2d_array(size);
  int ** vint=allocate_2d_array(size);
  int **interface;

  evolve_table_pm_boundary(size,table,beta, niter);

  for (long i=0; i<measurements; i++){
    interface=calc_interface(size, table,hint, vint);  
    add_crossing_number(size,hint,res);
    free_2d_array(size,interface);
    evolve_table_pm_boundary(size,table,beta,50);
  }
  free_2d_array(size,table);
  free_2d_array(size,hint);
  free_2d_array(size,vint);
  //  return res;
}

double measure_magnetisation(int size, double beta, int measurements)
{
  long long niter=3*size*size;
  int **table=allocate_2d_array(size);
  init_table_no_boundary(size,table);
  evolve_table_cyclic_boundary(size,table,beta, niter);
  double M=0;
  double mc=0;
  for (long i=0; i<measurements; i++){
    mc=calc_magnetization(size,table);
    M=M+mc/(size*size);
    evolve_table_cyclic_boundary(size,table,beta,50);
  }

  char name[100];
  sprintf(name,"d-%lf.xpm",1/beta);
  FILE *F=fopen(name,"w");
  print_table_xpm_simple(F,size,table);
  fclose(F);

  return M/(measurements);
}

double measure_susceptibility(int size, int **table, double beta, int measurements)
{
  long long niter=10000; //size*size*floor(1/fabs(1/beta-2.2698));

  evolve_table_cyclic_boundary(size,table,beta, niter);
  double M=0;
  double Msq=0;
  double mc=0;
  double Mfourth=0;  
  for (long i=0; i<measurements; i++){
    mc=((double)calc_magnetization(size,table))/(size*size);
    M=M+fabs(mc);
    Msq=Msq+mc*mc;
    Mfourth=Mfourth+mc*mc*mc*mc;

    evolve_table_cyclic_boundary(size,table,beta,1);
    printf("%lf ",mc);
    
  }
  M=M/measurements;
  Msq=Msq/measurements;
  Mfourth=Mfourth/measurements;
  
  free_2d_array(size,table);
  //  printf("%lf %lf %lf %lf %lf \n", 1/beta, M, beta*(Msq-M*M), beta*Msq,1.0-1.0/3.0*Mfourth/(Msq*Msq));
  return beta*(Msq-M*M);
}


double measure_susceptibility_cluster(int size, int **table, double beta, int measurements)
{
  long long niter=1000;
  long i;

  for (i=0;i<niter;i++) modify_cluster(size,table,rand_int(size), rand_int(size), beta);
  
  double M=0;
  double Msq=0;
  double mc=0;
  double Mfourth=0;  
  for (long i=0; i<measurements; i++){
    mc=((double)calc_magnetization(size,table))/(size*size);
    M=M+fabs(mc);
    Msq=Msq+mc*mc;
    Mfourth=Mfourth+mc*mc*mc*mc;

    modify_cluster(size,table,rand_int(size), rand_int(size), beta);

    printf("%lf ",mc);
    
  }
  M=M/measurements;
  Msq=Msq/measurements;
  Mfourth=Mfourth/measurements;
  
  free_2d_array(size,table);
  //  printf("%lf %lf %lf %lf %lf \n", 1/beta, M, beta*(Msq-M*M), beta*Msq,1.0-1.0/3.0*Mfourth/(Msq*Msq));
  return beta*(Msq-M*M);
}

double measure_Ul(int size,  double beta, int measurements)
{
  long long niter=size*size*size;
  int **table=allocate_2d_array(size);
  init_table_no_boundary(size,table);

  evolve_table_cyclic_boundary(size,table,beta, niter);
  double M=0;
  double Msq=0;
  double Mfourth=0;
  double mc=0;
  for (long i=0; i<measurements; i++){
    mc=((double)calc_magnetization(size,table))/(size*size);
    M=M+mc;
    Msq=Msq+mc*mc;
    Mfourth=Mfourth+mc*mc*mc*mc;
    evolve_table_cyclic_boundary(size,table,beta,1);
  }
  free_2d_array(size,table);
  M=M/measurements;
  Msq=Msq/measurements; 
  Mfourth=Mfourth/measurements;

  double res=1.0-1.0/3.0*Mfourth/(Msq*Msq);
  return res;
}

void print_Ul(double T, double step, long steps, long measures, long size1, long size2) {
  printf("U%ld_%ld:[",size1,size2);
  double Ul1,Ul2;
  Ul1=measure_Ul(size1,1/T,measures);
  Ul2=measure_Ul(size2,1/T,measures);
  for (int i=0;i<steps;i++) {
    printf("[%lf,%lf], ", T,Ul1/Ul2);
    T=T+step;
    Ul1=measure_Ul(size1,1/T,measures);
    Ul2=measure_Ul(size2,1/T,measures);
  }
  printf("[%lf,%lf] ", T,Ul1/Ul2);
  printf("]\n");
}

int** init_table_with_p(int size, double p){
  int **table=allocate_2d_array(size);
  long num=(long)((1.0-p)*size*size);
  long i=0,x,y;
  init_table_no_boundary(size,table);
  while (i<num) {
    x=rand_int(size);
    y=rand_int(size);
    if (table[x][y]!=0) {
      table[x][y]=0;
      i++;
    }
  }
  //  printf("%ld %ld %lf %lf\n", (long)size*size,num, (double)num/(size*size),p);
  return table;
}

int main(int argc, char **argv)
{



  gmp_randinit_mt(state);
  gmp_randseed_ui(state, time(NULL));
  mpf_init2(rop,256);
  
  double T=1.5;
  double step=0.01;
  int size=50;
  int steps = 100;
  long measures=10000;
  double p=1.0;
  char *method="";
  if (argc>1) size= atoi(argv[1]);
  if (argc>2) T=atof(argv[2]);
  if (argc>3) step=atof(argv[3]);
  if (argc>4) steps=atoi(argv[4]);
  if (argc>5) measures=atoi(argv[5]);
  if (argc>6) p=atof(argv[6]);
  if (argc>7) method=argv[7];


  //  printf("M:[\n");
  double M=0;
  if ((strlen(method)==4) && (strcmp("wolf",method)==0))
    for (int i=0;i<steps;i++) {
      T=T+step;
      printf("%lf ",T);
      measure_susceptibility_cluster(size, init_table_with_p(size,p),1/T,measures);
      printf("\n");
    }
  else
    for (int i=0;i<steps;i++) {
      T=T+step;
      printf("%lf ",T);
      measure_susceptibility(size, init_table_with_p(size,p),1/T,measures);
      printf("\n");
    }
  //  printf("[%lf,%lf] ", T,M);
  //  printf("]\n");
  
  //  print_Ul(T,step,steps, measures, size,4*size);
  //  print_Ul(T,step,steps, measures, 2*size,3*size);
  mpf_clear(rop);
  return 0;
}


//  int main()
//  {
//    double T=1.0;
//    double step=0.05;
//    int size=50;
//    srand(time(NULL));
//  //  FILE *F=fopen("debug.xpm","w");
//  //  int **table=allocate_2d_array(size);
//  //  evolve_table_cyclic_boundary(size,table,1/T,20*size*size);
//  //  print_table_xpm_simple(F,size,table);
//  //  free_2d_array(size,table);
//  //  fclose(F);
//    printf("[");
//    for (int i=0;i<80;i++) {
//      printf("[%lf,%lf], ", T,fabs(measure_magnetisation(size,1/T,1)));
//      T=T+step;
//    }
//    printf("]\n");
//    return 0;
//  }
//    

//int main()
//{
//  double beta=-0.440686854531;
//  int size=100;
//  long long niter=size*size*size;
//  int **table=allocate_2d_array(size);
//  int ** res = allocate_2d_array(size);
//  int ** hint=allocate_2d_array(size);
//  int ** vint=allocate_2d_array(size);
//  int **interface;
//
//  evolve_table(size,table,beta, niter);
//
//  interface=calc_interface(size, table,hint, vint);  
//
//  FILE *f;
//  f=fopen("table.xpm","w");
//  print_table_xpm_simple(f,size,table);
//  fclose(f);
//  f=fopen("hint.xpm","w");
//  print_table_xpm_simple(f,size,hint);
//  fclose(f);
//  f=fopen("vint.xpm","w");
//  print_table_xpm_simple(f,size,vint);
//  fclose(f);
//
//  add_crossing_number(size,vint,res);
//  f=fopen("res.xpm","w");
//  print_table_xpm_simple(f,size,res);
//  fclose(f);
//
//  free_2d_array(size,interface);
//  //  evolve_table(size,table,beta,10*size*size);
//
//  free_2d_array(size,table);
//  free_2d_array(size,hint);
//  free_2d_array(size,vint);
//  return 0;
//}
//  
//  int main(int argc, char **argv) 
//  {
//    int size=512;
//    double beta=0.440686854531;
//    long long niter=size*size*size;
//    FILE* Fout=stdout;
//  
//    if (argc>1) size= atoi(argv[1]);
//    if (argc>2) beta=atof(argv[2]);
//    if (argc>3) niter=atoll(argv[3]);
//    if (argc>4) Fout = fopen(argv[4],"w");
//  
//    srand(time(NULL));
//    
//  //  int **table=allocate_2d_array(size);
//  //
//  //  evolve_table(size,table,beta, niter);
//  //  //  int **interface = allocate_2d_array(size);
//  //  int **interface=calc_interface(size, table);
//  //
//  
//    long measurements=100;
//    int ** table= allocate_2d_array(size);
//    for (int i=0; i<20; i++)
//      compute_probability(size,beta,measurements,table);
//    print_table_as_list(Fout, size, 20*measurements, table);
//    fclose(Fout);
//  
//    free_2d_array(size,table);
//    //  free_2d_array(size,interface);
//    return 0;
//  }
