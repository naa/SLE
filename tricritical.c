#include "common.h"

int** init_table_with_delta(int M,int N)
{
  int **table=allocate_2d_rectangle(M,N);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) { 
      table[i][j] = rand_int(3)-1;
      //table[i][j] = 1;//-1+2*((i+j)%2);
    }
  }
  return table;
}

int modify_cell (int M, int N,int **table, int i, int j, double beta, double delta)
{
  int k,ns=0;
  int s_old=table[i][j];
  for (k=0; k<lattice_neighbors; k++)
    ns+=table[mod(i+lattice_dx[k],M)][mod(j+lattice_dy[k],N)];

  int s_new=(s_old+1+(nrand()<0.5));
  if (s_new>1) s_new=s_new-3;
  int de=-s_new*ns+delta*s_new*s_new-(-s_old*ns+delta*s_old*s_old);
  if ((de<0) || (nrand()<exp(-beta*de)))
    table[mod(i,M)][mod(j,N)]=s_new;
  return 0;
}


//int modify_cell (int size,int **table, int i, int j, double beta, double delta)
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
void modify_cluster (int M, int N,int **table, int i, int j, double beta)
{
  if (table[mod(i,M)][mod(j,N)]==0) return;   
  int *qx=(int*)malloc(M * N * sizeof(int));
  int *qy=(int*)malloc(M * N * sizeof(int));  
  long head=0;					
  long tail=0;
  int x,y,nx,ny;
  long blength=0;
  int k;
  int cluster=table[mod(i,M)][mod(j,N)];
  table[mod(i,M)][mod(j,N)]=3;
  qx[tail]=i;
  qy[tail]=j;
  tail++;
  while (tail>head) {
    x=qx[head];
    y=qy[head];
    head++;
    for (k=0; k<lattice_neighbors; k++) {
      nx=mod(x+lattice_dx[k],M);
      ny=mod(y+lattice_dy[k],N);      
      if ((table[nx][ny]==cluster) && (nrand()<1.0-exp(-2*beta))){
	qx[tail]=nx;
	qy[tail]=ny;
	table[nx][ny]=3;
	tail++;
      } else if (table[nx][ny]==-cluster) {
	blength++;
      }
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

void evolve_table_cyclic_boundary ( int M, int N, int** table, double beta, double delta, long long niter ) 
{
  int i,j;

  for (long long k = 0; k < niter; k++){
    for (i=0;i<M;i++)
      for (j=0;j<N;j++)
	//modify_cluster (size,table,rand_int(size-1), rand_int(size-1),beta);
	modify_cell (M, N,table,rand_int(M), rand_int(N),beta,delta);
	//modify_cell (size,table,i,j,beta);
  }
}

//void evolve_table_pm_boundary ( int M, int N, int** table, double beta, long long niter ) 
//{
//  int i,j;
//  init_table_pm_boundary(M,N,table);
//
//  for (long long k = 0; k < niter; k++){
//    i=1 + rand_int(M-2);
//    j=1 + rand_int(N-2);
//    modify_cell (M,N,table,i, j,beta);
//  }
//}
//
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
  fprintf(F, "\"%d %d %d %d\",\n", size, size, 3, 1);
  fprintf(F, "\"+ c #ffffff\",\n");
  fprintf(F, "\"0 c #ff0000\",\n");  
  fprintf(F, "\"- c #000000\",\n");

  for (int i = 0; i < size; i++) {
    fprintf(F, "\"");
    for (int j = 0; j < size; j++)
      if (table[i][j] > 0) 
	fprintf(F, "+"); 
      else if (table[i][j]<0)
	fprintf(F, "-");
      else
	fprintf(F, "0");	
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

//void compute_probability(int size, double beta, long measurements,int **res) 
//{
//  long long niter=3*size*size;
//  int **table=allocate_2d_array(size);
//  //  int ** res = allocate_2d_array(size);
//  int ** hint=allocate_2d_array(size);
//  int ** vint=allocate_2d_array(size);
//  int **interface;
//
//  evolve_table_pm_boundary(size,size,table,beta, niter);
//
//  for (long i=0; i<measurements; i++){
//    interface=calc_interface(size, table,hint, vint);  
//    add_crossing_number(size,hint,res);
//    free_2d_array(size,interface);
//    evolve_table_pm_boundary(size,size,table,beta,50);
//  }
//  free_2d_array(size,table);
//  free_2d_array(size,hint);
//  free_2d_array(size,vint);
//  //  return res;
//}
//
//double measure_magnetisation(int size, double beta, int measurements)
//{
//  long long niter=3*size*size;
//  int **table=allocate_2d_array(size);
//  init_table_no_boundary(size,size,table);
//  evolve_table_cyclic_boundary(size,size,table,beta, niter);
//  double M=0;
//  double mc=0;
//  for (long i=0; i<measurements; i++){
//    mc=calc_magnetization(size,size,table);
//    M=M+mc/(size*size);
//    evolve_table_cyclic_boundary(size,size,table,beta,50);
//  }
//
//  char name[100];
//  sprintf(name,"d-%lf.xpm",1/beta);
//  FILE *F=fopen(name,"w");
//  print_table_xpm_simple(F,size,table);
//  fclose(F);
//
//  return M/(measurements);
//}
//

double calc_energy(int M,int N,int **table,double delta) {
  int i, j,k;
  double e=0;

  for (i=0; i<M; i++) 
    for (j=0;j<N; j++) {
        for (k=0; k<lattice_neighbors; k++)
	  e-=table[i][j]*table[mod(i+lattice_dx[k],M)][mod(j+lattice_dy[k],N)];
	e+=2*delta*table[i][j]*table[i][j];
    }
  return e/2;
}
  

double measure_susceptibility(int M, int N, int **table, double beta, double delta, int measurements)
{
  //  long long niter=10000; //size*size*floor(1/fabs(1/beta-2.2698));
  long long niter=10000;
  evolve_table_cyclic_boundary(M,N,table,beta,delta, niter);
  double mm=0;
  double Msq=0;
  double mc=0;
  double Mfourth=0;
  double E=0;
  double EE=0;
  double Esq=0;
  double cv=0;
  fprintf(stderr,"Thermalization completed: T=%lf, d=%lf\n",1/beta,delta);  
  for (long i=0; i<measurements; i++){
    mc=((double)calc_magnetization(M,N,table))/(M*N);
    mm=mm+fabs(mc);
    Msq=Msq+mc*mc;
    Mfourth=Mfourth+mc*mc*mc*mc;
    E=calc_energy(M,N,table,delta)/(M*N);
    //    fprintf(stderr,"%lf\n",E);
    EE=EE+E;
    Esq=Esq+E*E;
    evolve_table_cyclic_boundary(M,N,table,beta,delta,1);
    //    printf("%lf ",mc);
    
  }
  mm=mm/measurements;
  Msq=Msq/measurements;
  Mfourth=Mfourth/measurements;
  EE=EE/measurements;
  Esq=Esq/measurements;
  cv=beta*beta*(Esq-EE*EE);
  free_2d_array(M,table);
  printf("%lf %lf %lf %lf %lf %lf %lf \n", 1/beta, mm, beta*(Msq-mm*mm), beta*Msq,1.0-1.0/3.0*Mfourth/(Msq*Msq),EE,cv);
  return beta*(Msq-mm*mm);
}


double measure_susceptibility_cluster(int M, int N, int **table, double beta, int measurements)
{
  long long niter=1000;
  long i;

  for (i=0;i<niter;i++) modify_cluster(M,N,table,rand_int(M), rand_int(N), beta);
  
  double mm=0;
  double Msq=0;
  double mc=0;
  double Mfourth=0;
  long E=0;
  for (long i=0; i<measurements; i++){
    mc=((double)calc_magnetization(M,N,table))/(M*N);
    mm=mm+fabs(mc);
    Msq=Msq+mc*mc;
    Mfourth=Mfourth+mc*mc*mc*mc;

    modify_cluster(M,N,table,rand_int(M), rand_int(N), beta);
    //    E=state_energy(M,N,table);
    //    printf("%lf %ld\n",mc, E);
    
  }
  mm=mm/measurements;
  Msq=Msq/measurements;
  Mfourth=Mfourth/measurements;
  
  free_2d_array(M,table);
  printf("%lf %lf %lf %lf %lf \n", 1/beta, mm, beta*(Msq-mm*mm), beta*Msq,1.0-1.0/3.0*Mfourth/(Msq*Msq));
  return beta*(Msq-mm*mm);
}

//double measure_Ul(int size,  double beta, int measurements)
//{
//  long long niter=size*size*size;
//  int **table=allocate_2d_array(size);
//  init_table_no_boundary(size,size,table);
//
//  evolve_table_cyclic_boundary(size,size,table,beta, niter);
//  double M=0;
//  double Msq=0;
//  double Mfourth=0;
//  double mc=0;
//  for (long i=0; i<measurements; i++){
//    mc=((double)calc_magnetization(size,size,table))/(size*size);
//    M=M+mc;
//    Msq=Msq+mc*mc;
//    Mfourth=Mfourth+mc*mc*mc*mc;
//    evolve_table_cyclic_boundary(size,size,table,beta,1);
//  }
//  free_2d_array(size,table);
//  M=M/measurements;
//  Msq=Msq/measurements; 
//  Mfourth=Mfourth/measurements;
//
//  double res=1.0-1.0/3.0*Mfourth/(Msq*Msq);
//  return res;
//}

//void print_Ul(double T, double step, long steps, long measures, long size1, long size2) {
//  printf("U%ld_%ld:[",size1,size2);
//  double Ul1,Ul2;
//  Ul1=measure_Ul(size1,1/T,measures);
//  Ul2=measure_Ul(size2,1/T,measures);
//  for (int i=0;i<steps;i++) {
//    printf("[%lf,%lf], ", T,Ul1/Ul2);
//    T=T+step;
//    Ul1=measure_Ul(size1,1/T,measures);
//    Ul2=measure_Ul(size2,1/T,measures);
//  }
//  printf("[%lf,%lf] ", T,Ul1/Ul2);
//  printf("]\n");
//}
//


int main(int argc, char **argv)
{
  lattice_neighbors=square_lattice_neighbors;
  lattice_dx=square_lattice_dx;
  lattice_dy=square_lattice_dy;
  
  gmp_randinit_mt(state);
  gmp_randseed_ui(state, time(NULL));
  mpf_init2(rop,256);
  //  wang_landau(5,25,1.0,20,0.8,5,1.5, 0.05, 30);
  //  return 0;

//  FILE *f;
//  f=fopen("tric-table.xpm","w");
//  int ** tab=init_table_with_delta(20,20);
//  print_table_xpm_simple(f,20,tab);
//  evolve_table_cyclic_boundary(20,20,tab,1/2.25,0.0, 100);
//  print_table_xpm_simple(f,20,tab);  
//  fprintf(f,"\n%lf %lf\n",(double)calc_magnetization(20,20,tab),calc_energy(20,20,tab,0.0));
//  fclose(f);
  
  double T=1.5;
  double step=0.01;
  int M=50;
  int N=50;
  int steps = 100;
  long measures=10000;
  double delta=1.0;
  char *method="";

  int flags, opt;
  int nsecs, tfnd;

  while ((opt = getopt(argc, argv, "M:N:T:S:s:d:m:a:L:")) != -1) {
    switch (opt) {
    case 'M':
      M=atoi(optarg);
      break;
    case 'N':
      N=atoi(optarg);
      break;
    case 'T':
      T = atof(optarg);
      break;
    case 'S':
      step = atof(optarg);
      break;
    case 's':
      steps = atoi(optarg);
      break;
    case 'd':
      delta = atof(optarg);
      break;
    case 'm':
      measures=atol(optarg);
      break;
    case 'a':
      method=optarg;
      break;
    case 'L':
      if ((strlen(optarg)==6) && (strcmp("triang",optarg)==0)) {
	lattice_neighbors=triangular_lattice_neighbors;
	lattice_dx=triangular_lattice_dx;
	lattice_dy=triangular_lattice_dy;
      }
      break;
      
    default: /* '?' */
      fprintf(stderr, "Usage: %s [-M xsize] [-N ysize] [-T temperature] [-S temperature step] [-s number of steps] [-d delta] [-m measures] [-a algorithm:metropolis|wolf] [-L lattice:square|triang]  name\n",
	      argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  

  //  printf("M:[\n");
  if ((strlen(method)==4) && (strcmp("wolf",method)==0))
    for (int i=0;i<steps;i++) {
      T=T+step;
      //      printf("Temperature: %lf \n",T);
      measure_susceptibility_cluster(M,N, init_table_with_delta(M,N),1/T,measures);
      //      printf("\n");
    }
  else
    for (int i=0;i<steps;i++) {
      T=T+step;
      measure_susceptibility(M,N, init_table_with_delta(M,N),1/T,delta,measures);
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
