#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){
  int *c, d;
  d = 10;
  c = &d;
  printf("c is : %d\n", c);
  printf("d is : %d\n", d);
  d = *c + 1;
  printf("d is : %d\n", d);

}
