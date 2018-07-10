#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int combinations ( int n, int k, unsigned long long int *b,int *tab);
unsigned long long int choose(int n, int k);

int main(){
  int *table, sites, electrons, dimension, i;
  unsigned long long int *b;
  sites = 9;
  electrons = 2;
  dimension=choose(sites,electrons);
  printf("dim = %d\n", dimension);
  //table=(int*) malloc(electrons*dimension*sizeof(int));
  printf("malloc input = %d\n", (electrons*dimension*sizeof(int)));

  table= malloc(electrons*dimension*sizeof(int));
  b = malloc(dimension);
  printf("3left is : %d\n", (1ULL<<10));
  combinations (sites,electrons,b,table);
  for (i=0;i<dimension;i++) printf("p    = %p\n", b[i]);
  printf("TABLEINCOMING\n");
  for (i=0;i<(dimension+1);i++) printf("p    = %p\n", table[i]);


  //printf("p    = %p\n", (void *) b);

}

unsigned long long int choose(int n, int k) //the number of combination of k out of n
{
     int i;
     unsigned long long int c;
     c=1ULL;
     for (i=0;i<k;i++) c=c*(n-i);
     for (i=0;i<k;i++) c=c/(i+1);
     return c;
}

int combinations ( int n, int k, unsigned long long int *b,int *tab)
{
     unsigned long long int x,y;
     int i,c,d;
     x=0ULL;
     for (i=0;i<k;i++)
     {
         x=x+(1ULL<<i);
     }
     b[0]=x;
     c=0;
     d=0;
     i=0;
     while ((c<n)&& (d<k))
     {
         if (x & (1ULL<<c))
         {
             tab[i*k+d]=c+1;
             d++;
         }
         c++;
     }
     for (i=1;i<choose(n,k);i++)
     {
         y = (x | (x - 1)) + 1;
         x = y | ((((y & -y) / (x & -x)) >> 1) - 1);
         b[i]=x;
         c=0;
         d=0;
         while ((c<n)&& (d<k))
         {
             if (x & (1ULL<<c))
             {
                 tab[i*k+d]=c+1;
                 d++;
             }
             c++;
         }
     }
}
