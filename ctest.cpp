#include <stdio.h>
#include <stdlib.h>
#include <string.h>


dostuff(){
  int bit, j;
  bit = 1ULL;
  j = (bit << 2) + (bit <<4);
  printf("%i", j);
}

int main(){
  dostuff();
}
