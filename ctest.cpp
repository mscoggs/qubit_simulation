#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define nx 4


int main(){
  printf(   "%i"     ,!((1ULL<<2)&4));
  printf(   "%i"     ,!((1ULL<<3)&4));
  printf(   "%i"     ,!((1ULL<<1)&4));
  printf(   "%i"     ,((1ULL<<1)&4));
  printf(   "%i\n"     ,((1ULL<<2)&4));
  printf(   "%i"     ,(1ULL<<6));
  if(1) printf("True");
  if(0) printf("construct_model_hamiltonian");
  if(7) printf("Truen");
}
