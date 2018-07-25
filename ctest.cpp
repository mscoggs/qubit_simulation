#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define nx 4

void construct_lattice(int lattice[][nx]){
  int x,i;
  for (x=0;x<nx;x++)
  {
    if (x%2 ==1)
    {
      for (i=0;i<nx;i++)
      {
        lattice[x][i] = (nx*x)+i+1;
      }
    }
    else
    {
      for (i=0;i<nx;i++)
      {
        lattice[x][nx-1-i] = nx*x+i+1;
      }
    }
  }
}

int main(){
  int i,x;
  //int* lattice[nx][nx];
  //lattice=(int*) malloc(nx*nx*sizeof(int));
  int lattice[nx][nx];
  construct_lattice(lattice);
  for (i=0;i<nx;i++)
  {
    for (x=0;x<nx;x++)
    {
      printf("%i\n", lattice[i][x]);
    }
  }
//  for (x=0;x<nx;x++)

//  for
}
