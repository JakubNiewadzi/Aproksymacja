#include "makespl.h"
#include "gaus/piv_ge_solver.h"
#include "gauss.h"


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>



double makeHermite (double x, int n)  //tworzy funkcje do bazy
{    
  if(n==0){
    return 1;
  }else if(n==1){
    return 2*x;
  }else if(n<0){
    return 0;
  }else{
  double n1=1;
  double n2=2*x;
  double n3;
  for(int i=1;i<n;i++){
    n3=2*x*n2-2*i*n1;
    n1=n2;
    n2=n3;
  }
  return n2;
  }
}

void multiplication(int m, int o, int n, matrix_nasz * tab1, matrix_nasz * tab2, matrix_nasz * nowy){
  double ile;
  for (int i=0;i<m;i++)   // m to wysokość(ilość wierszy) pierwszej tablicy
    {
    for(int j=0;j<o;j++)  // o to szerokość(ilość kolumn) drugiej tablicy
    {
      ile=0;
      for(int k=0;k<n;k++)  // n to wysokość drugiej i szerokość pierwszej tablicy
    {
        ile+=tab1->e[i][k]*tab2->e[k][j];  
    }
      nowy->e[i][j]=ile; //tab3 to tablica z wynikiem
  }
}
}


void our_make(points_t* pts, char* out, double fromX, double toX, int n) {
  //matrix_nasz D=make_matrix(,x = pts->x);
  FILE *fp;
  if(out !=NULL){
    fp = fopen(out, "w");
  }else{
    fp = fopen("wyjscie.txt", "w");
  }
  matrix_nasz *D= NULL;
  matrix_nasz *Dt;
  double *A;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
  int		i, j, k;
	int		nb = pts->n - 3 > 9 ? 9 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	//D = make_matrix(pts->n, nb); ------------------------------------------------
  D=malloc(sizeof(*D));
  D->e=malloc(pts->n*sizeof(double*));
  for(int i=0;i<pts->n;i++){
    D->e[i]=malloc(nb*sizeof(double));
  }
  for(int i=0;i<pts->n;i++)
  {
      for(int j=0;j<nb;j++)
    {
      D->e[i][j]=makeHermite(x[i],j);
      //printf("%lf ",D->e[i][j] );
    }
    //printf("\n");
  }
  //printf("HERMIT TU %lf \n",makeHermite(150,2));
  
  A=malloc(nb*sizeof(double));
  matrix_nasz * F=malloc(sizeof (*F));
  Dt = malloc(sizeof(*Dt));
  Dt->e=malloc(nb*sizeof(double*));
  for(int i=0;i<nb;i++){
    Dt->e[i]=malloc(pts->n*sizeof(double));
  }
  for(int i=0;i<nb;i++)
  {
      for(int j=0;j<pts->n;j++)
    {
      Dt->e[i][j]=D->e[j][i];
      //printf("%lf ",Dt->e[i][j] );
    }
    //printf("\n");
  }
  //printf("\n\n");
  matrix_nasz * DtD=malloc(sizeof *DtD);
  DtD->e=malloc(nb*sizeof(double*));
  for(int i=0;i<nb;i++){
    DtD->e[i]=malloc(nb*sizeof(double));
  }
  multiplication(nb, nb, pts->n, Dt, D , DtD);
  for(int i=0;i<nb;i++){
    for(int j=0;j<nb;j++){
      //printf("%.2lf ", DtD->e[i][j]);
    }
    //printf("\n");
  }
  //printf("pts %d",pts->n);

  F->e=malloc(pts->n*sizeof(double*));
  for(int i=0;i<pts->n;i++){
    F->e[i]=malloc(sizeof(double));
    F->e[i][0]=pts->y[i];
  }
  
  matrix_nasz * DtF=malloc(sizeof (*DtF));
  DtF->e=malloc(nb*sizeof(double*));
  for(int i=0;i<nb;i++){
    DtF->e[i]=malloc(sizeof(double));
  }
  multiplication(nb, 1,  pts->n, Dt, F ,DtF);

  for(int i=0;i<nb;i++){
    for(int j=0;j<1;j++){
     // printf("%.2lf ", DtF->e[i][j]);
    }
   // printf("\n");
  }
  double *wsp = malloc(nb*sizeof(double));
  wsp=gauss(nb, DtD, DtF, wsp);

  //for (int i=0;i<nb;i++)
    //printf("%lf ",wsp[i]);
  //printf("\n");
  double wart=0;
  //printf("wartosci1:\n");
  for (int j=0;j<pts->n;j++) //dla każdego x'a
  {
    for (int i=0;i<nb;i++)  //dla każdej wartości hermita
    {
      //printf("%lf razy %lf", wsp[i], makeHermite(pts->x[j], i));
    wart+=wsp[nb-i-1]*makeHermite(pts->x[j],i);
    }
    printf("dla x=%lf",pts->x[j]);
    printf(" dostajemy %lf a powinnismy %lf \n",wart,pts->y[j]);
    
    wart=0;
  }
  double start;
  double koniec;
  double skok;
  int max;
  if(fromX!=0){
    start=fromX;
  }else
  {
    start=pts->x[0];
  }
  if(toX!=0){
    koniec=toX;
  }else
  {
    koniec=pts->x[pts->n-1];
  }
  if(n>1){
    skok=(koniec-start)/(n-1);
    max=n;
  }else
  {
    skok=(koniec-start)/(10000-1);
    max = 10000;
  }
  printf("%s %lf, %lf, %d, %lf", out, start, koniec, n, skok);
  for(int j=0;j<max;j++)
  {
    for (int i=0;i<nb+2;i++)  //dla każdej wartości hermita
    {
      //printf("%lf razy %lf", wsp[i], makeHermite(pts->x[j], i));
    wart+=wsp[nb-i-1]*makeHermite(start+skok*j,i);
    }
   fprintf( fp, "%lf %lf\n", start+skok*j, wart);
   wart=0;
  }
  for (int j=0;j<pts->n;j++) //dla każdego x'a
    {
      for (int i=0;i<nb;i++)  //dla każdej wartości hermita
      {
        //printf("%lf razy %lf", wsp[i], makeHermite(pts->x[j], i));
      wart+=wsp[nb-i-1]*makeHermite(pts->x[j],i);
      }
      //printf("dla x=%lf",pts->x[j]);
      //printf(" dostajemy %lf a powinnismy %lf \n",wart,pts->y[j]);
      
      wart=0;
    }
    
  fclose(fp);
  /*
  for(int i=0;i<pts->n;i++){ 
    free(D->e[i]);
  }
  free(D->e);
  free(D);
  
  for(int i=0;i<nb;i++){
    free(DtD->e[i]);
  }
  free(DtD->e);
  free(DtD);
  

  for(int i=0;i<nb;i++){
    free(DtF->e[i]);
  }
  free(DtF->e);
  free(DtF);
  
  
  free(wsp);
  

  for(int i=0;i<pts->n;i++){
    free(F->e[i]);
  }
  free(F->e);
  free(F);

  free(A);                

  for(int i=0;i<nb;i++){
    free(Dt->e[i]);
  }
 
  free(Dt->e);
  free(Dt);
  */
  }