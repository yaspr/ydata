#include <stdio.h>

#include "types.h"
#include "ydata.h"

i32 main(void)
{
  u64 n = 100;
  ydata_t *d1 = ydata_create("test", n);
  ydata_t *d2 = ydata_create(NULL, n);
    
  for (u64 i = 0; i < d1->n_datapoints; i++)
    {
      d1->datapoints[i] = (f64)(i + 2) / (f64)(i + 3);
      d2->datapoints[i] = (f64)(i + 4) / (f64)(i + 2);
    }

  ydata_analyze(d1);
  ydata_analyze(d2);
  
  f64 p = ydata_pearson_correlation(d1, d2);
  
  printf("Testing:\n\n");
  
  ydata_print(d1);
  ydata_print(d2);
  
  ydata_dump(d1, "data1.csv");
  ydata_dump(d2, "data2.csv");
  
  printf("\nPearson correlation factor: %lf\n", p);
  
  ydata_destroy(&d1);
  ydata_destroy(&d2);
}
