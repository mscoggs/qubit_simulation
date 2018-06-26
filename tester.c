unsigned long long int C(int n, int k) //the number of combination of k
out of n
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


     for (i=1;i<C(n,k);i++)
     {
         //y= x | (x - 1);
         //x = (y + 1) | (((~y & -~y) - 1) >> (__builtin_ctz(x) + 1));


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
