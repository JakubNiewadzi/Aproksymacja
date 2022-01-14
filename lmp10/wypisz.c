#include <stdio.h>


void wypisz(double **tab, double *wynik, int x){
   for (int j=0;j<x;j++)
    {
        for (int k=0;k<x;k++)
        {
            printf("%7.4lf ",tab[j][k]);
        }
        printf("%7.4lf \n",wynik[j]);
    }
    printf("\n");
}