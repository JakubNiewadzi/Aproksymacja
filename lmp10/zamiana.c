#include <stdio.h>
#include <math.h>


void zamien(double **tab, double *wynik, int i, int x, int *zmienione){
  double max = tab[i][0];
  double tmp;
  int pomoc;
  int max_j=0;
  for(int j=0;j<x-i;j++){
    if(fabs(tab[j][i])>fabs(max)){
      max=tab[j][i];
      max_j=j;
    }
  }
  pomoc=zmienione[max_j];
  zmienione[max_j]=zmienione[x-i-1];
  zmienione[x-i-1]=pomoc;
  
  if(max!=tab[i][x-1-i]){
    for(int j =0;j<x;j++){ 
      tmp = tab[x-1-i][j];
      tab[x-1-i][j]=tab[max_j][j];
      tab[max_j][j]=tmp;
    }
    tmp=wynik[x-1-i];
    wynik[x-1-i]=wynik[max_j];
    wynik[max_j]=tmp;
  }
}