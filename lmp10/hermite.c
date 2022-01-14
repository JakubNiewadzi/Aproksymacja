#include "makespl.h"
#include "gaus/piv_ge_solver.h"
#include "gauss.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>



double makeHermite(double x, int n){     //tworzy funkcje do bazy
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

matrix_t * multiplication(int m, int o, int n, matrix_t * tab1, matrix_t * tab2){
  int ile;
  matrix_t * tab3=make_matrix(n, n);

  tab3->e=malloc(n*sizeof(double*));
  for(int i=0;i<n;i++){
    tab3->e[i]=malloc(n*sizeof(double));
  }


  for (int i=0;i<m;i++)   // m to wysokość(ilość wierszy) pierwszej tablicy
    {
    for(int j=0;j<o;j++)  // o to szerokość(ilość kolumn) drugiej tablicy
    {
      ile=0;
      for(int k=0;k<n;k++)  // n to wysokość drugiej i szerokość pierwszej tablicy
    {
        ile+=tab1->e[i][k]*tab2->e[k][j];  
    }
      tab3->e[i][j]=ile; //tab3 to tablica z wynikiem
  }
}
  return tab3;
}






/*double f(hermite hermite, double x){  //
  double k=0;
  for(int i=0;i<hermite.nwsp;i++){
    k+=hermite.wsp[i]*makeHermite(x, i);
  }
  printf("%lf\n", k);
  return k;
}

double df(hermite hermite, double x, int level){
  double k=0;
  int m=1;
  int n=1;
  for(int i=0;i<hermite.nwsp;i++){
    m=1;
    n=1;
    if(i<2){
      m=1;
      n=1;
    }else if (i-level<2){
      for(int j=1;j<i;j++){
        n*=j+1;
      }
      m=1;
    }else{
      for(int j=1;j<i;j++){
        n*=j+1;
      }
      for(int j=1;j<i-level;j++){
        m*=j+1;
      }
    }
    k+=2*(n/m)*hermite.wsp[i]*makeHermite(x, i-level);
  }
  printf("pochodna: %lf\n", k);
  return k;
}*/


void make_spl(points_t* pts, spline_t* spl) {
  //matrix_t D=make_matrix(,x = pts->x);
  
  matrix_t       *D= NULL;
  matrix_t *Dt;
  double *A;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
  int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	D = make_matrix(pts->n, nb);
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
    printf("\n");
  }
  
  
  A=malloc(nb*sizeof(double));
  matrix_t * F=make_matrix(pts->n, 1);
  Dt = make_matrix(nb,pts->n);
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
  matrix_t * DtD=make_matrix(nb, nb);
  DtD->e=malloc(nb*sizeof(double*));
  for(int i=0;i<nb;i++){
    DtD->e[i]=malloc(nb*sizeof(double));
  }
  DtD=multiplication(nb, nb, pts->n, Dt, D );
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
  
  matrix_t * DtF=make_matrix(nb, 1);
  DtF->e=malloc(nb*sizeof(double*));
  for(int i=0;i<nb;i++){
    DtF->e[i]=malloc(sizeof(double));
  }
  DtF=multiplication(nb, 1,  pts->n, Dt, F );

  for(int i=0;i<nb;i++){
    for(int j=0;j<1;j++){
      printf("%.2lf ", DtF->e[i][j]);
    }
    printf("\n");
  }
  double *wsp = malloc(nb*sizeof(double));
  wsp = gauss(nb, DtD, DtF);


  double wart=0;
  for (int j=0;j<pts->n;j++)
  {
    for (int i=0;i<nb;i++)
    {
    wart+=wsp[i]*makeHermite(pts->x[j],i);
    }
    printf("%lf \n",wart);
    wart=0;
  }
  // WAŻNE!!!!!!!!!!!!!
/*
  //Mnożenie macierzy
  //Chcemy to zrobić w funkcji czy oddzielnie dla każdego?
int ile;
for (int i=0;i<m;i++)   // m to wysokość(ilość wierszy) pierwszej tablicy
{
  for(int j=0;j<o;j++)  // o to szerokość(ilość kolumn) drugiej tablicy
  {
    ile=0;
    for(int k=0;k<n;k++)  // n to wysokość drugiej i szerokość pierwszej tablicy
    {
      ile+=tab1[i][k]*tab2[k][j]  
    }
    tab3[i][j]=ile; //tab3 to tablica z wynikiem
  }
}
*/
//Do zrobienia: Dt * D * A = Dt * f, w sumie 2 mnożenia i rozwiązanie Gaussa

  
  
  
  
/*
eqs = make_matrix(nb, nb + 1);
	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

  hermite hermite;
  hermite.nwsp=nb;  //liczba wspolczynnikow
  hermite.wsp=malloc(sizeof(double) * nb);
  hermite.a=a;      //poczatek zakresu
  hermite.b=b;      //koniec zakresu

  
  for (j = 0; j < nb; j++) {
	  for (i = 0; i < nb; i++)
		  for (k = 0; k < pts->n; k++)
			  add_to_entry_matrix(eqs, j, i, f(hermite, x[k])*f(hermite, x[k]));
        printf("%lf\n", eqs->e[k]);

	  for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k]*f(hermite, x[k]));
	}
  

  if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * f  (hermite, xx);
				spl->f1[i] += ck * df (hermite, xx, 1);
				spl->f2[i] += ck * df (hermite, xx, 2);
				spl->f3[i] += ck * df (hermite, xx, 3);
			}
		}
	}*/
  double a1=2;
  long double b1=2;
  for (int i=0;i<10;i++)
  {
    a1=a1*a1;
    b1=b1*b1;
    printf("%lf vs %Lf \n",a1,b1);
  }








}