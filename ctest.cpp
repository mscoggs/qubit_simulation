#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define nx 4


int main(){
  unsigned long long int p,i, z_mask,*a;
  unsigned long long int vtemp;
  int z_count = 0;
  if((1ULL<<4) & (1ULL<<4))
    z_count += 7;
  if((1ULL<<3) & (1ULL<<4))
    z_count += 5;
  printf("%i", z_count);

}
