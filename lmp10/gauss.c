#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#define BUFSIZE 8192

#include "gaus/piv_ge_solver.h"
#include "wypisz.h"
#include "zamiana.h"
#include "sposob.h"

double * gauss(int n,matrix_t * rownania, matrix_t * wartosci) //rownania to DtD a wartosci to DtF n to glebokosc
{
    int error=0;

    
   /* 
    FILE *in= argc > 1 ? fopen( argv[1], "r" ) : stdin;
    if(fscanf(in, "%d\n", &x)!=1){
      fprintf(stderr, "Ups! Podano błędny format danych!");
      return 1;
    }
  */


    double **tab;
    double *wynik;
    double *zmienne;
    int *zmienione;
    zmienne = malloc( n * sizeof( double ) );
    wynik = malloc( n * sizeof( double  ) );
    tab = malloc( n * sizeof( double * ) );
    zmienione = malloc( n * sizeof( int  ) );

    for (int i = 0; i < n; i++ )
    {
      zmienione[i]=i;
    }
  
    for (int i = 0; i < n; i++ )
    {
        tab[i] = malloc( n * sizeof( double ) );
    }
    
    tab=rownania->e;

    for(int i=0;i<n;i++){
      wynik[i]=wartosci->e[i][0];
    }


    // ZAMIANA -----------------------------------------

   //wypisz(tab, wynik, n);
  

  int pom=n-1;
  double przezco;
  for (int i = 0; i < n-1; i++ )
  {
    zamien(tab, wynik, i, n, zmienione);
    for (int j=0;j<n-1-i;j++) //zerowanie kolejnych wierszy
    {
      double przezco=tab[j][i]/tab[pom-i][i];
      for (int k=i;k<n;k++) // mnożenie wiersza przez wyznaczoną stałą
        {
          tab[j][k]-=(tab[pom-i][k]*przezco);        
          
           //POMNOŻYĆ WYNIK!!!
        }
      wynik[j]-=przezco*wynik[pom-i];
  
    }
    
  }
  

  for (int i = 0; i < n; i++ )
  {
    for (int j=0;j<i;j++)
    {
      wynik[i]-=zmienne[j]*tab[i][n-1-j];
    }
    zmienne[i]=wynik[i]/tab[i][n-i-1];
  }

  if(spr(n, tab)!=0){
    fprintf(stderr, "Ups! Podana macież jest osobliwa i nie mogę rozwiącać równania.");
  }

  printf("Wyniki:\n");
  for (int i = 0; i < n; i++ )
  {
    printf("x%d: %lf\n",zmienione[i]+1, zmienne[i]);
  }
  double* dzialaj=malloc(sizeof(double)*n);
  for (int i=0;i<n;i++)
  {
    dzialaj[zmienione[i]]=zmienne[i];
  }
  
  printf("\n");
for (int i = 0; i < n; i++ )
  {
    printf("x%d: %lf\n",i+1, dzialaj[i]);
  }

  return dzialaj;

}