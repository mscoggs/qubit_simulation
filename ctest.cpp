#include <stdio.h>
#include </usr/local/include/gsl/gsl_rng.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  int i,;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  double u = gsl_rng_uniform (r);
  printf ("%.5f\n", u);

  gsl_rng_free (r);

  return 0;
}
