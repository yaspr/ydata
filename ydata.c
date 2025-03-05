#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "ydata.h"
#include "../ynotif/ynotif.h"

void _sort(f64 *datapoints, u64 n)
{
  for (u64 i = 0; i < n; i++)
    for (u64 j = i + 1; j < n; j++)
      if (datapoints[i] > datapoints[j])
	{
	  f64 tmp = datapoints[i];
	  
	  datapoints[i] = datapoints[j];
	  datapoints[j] = tmp;
	}
}

f64 _mean(f64 *datapoints, u64 n)
{
  f64 s = 0.0;
  
  for (u64 i = 0; i < n; i++)
    s += datapoints[i];
  
  return s / (f64)n;
}

f64 _variance(f64 *datapoints, u64 n, f64 mean)
{
  f64 v = 0.0;
  
  for (u64 i = 0; i < n; i++)
    v += (datapoints[i] - mean) * (datapoints[i] - mean);
  
  v /= (f64)(n - 1);
  
  return v;
}

f64 _pearson(f64 *datapoints1, f64 *datapoints2, u64 n)
{
  f64 sum_a  = 0.0;
  f64 sum_b  = 0.0;
  f64 sum_a2 = 0.0;
  f64 sum_b2 = 0.0;
  f64 sum_ab = 0.0;
  f64 nn     = (f64)n;
  
  for (u64 i = 0; i < n; i++)
    {
      sum_a += datapoints1[i];
      sum_b += datapoints2[i];

      sum_a2 += (datapoints1[i] * datapoints1[i]);
      sum_b2 += (datapoints2[i] * datapoints2[i]);

      sum_ab += (datapoints1[i] * datapoints2[i]);
    }
  
  f64 top    = (nn * sum_ab) - (sum_a * sum_b);
  f64 bottom = sqrt((nn * sum_a2) - (sum_a * sum_a)) * sqrt((nn * sum_b2) - (sum_b * sum_b));
  
  return (top / bottom);
}

ydata_t *ydata_create(const ascii *title, u64 n_datapoints)
{
  if (n_datapoints < YDATA_MIN_DATAPOINTS)
    ynotif_error("ydata_create", YNOTIF_EXIT, "cannot allocate memory for data\n");
  
  ydata_t *d = malloc(sizeof(ydata_t));
  
  if (!d)
    ynotif_error("ydata_create", YNOTIF_EXIT, "cannot allocate memory for data\n");
  
  d->title = title;
  
  d->n_datapoints = (n_datapoints & 1) ? n_datapoints : n_datapoints + 1;
  
  d->datapoints = malloc(sizeof(f64) * d->n_datapoints);

  if (!d->datapoints)
    ynotif_error("ydata_create", YNOTIF_EXIT, "cannot allocate memory for datapoints\n");
  
  for (u64 i = 0; i < d->n_datapoints; i++)
    d->datapoints[i] = 0.0;
  
  d->min                      = 0.0;
  d->max                      = 0.0;
  d->mean                     = 0.0;
  d->median                   = 0.0;
  d->variance                 = 0.0;
  d->min_median_diff          = 0.0;
  d->min_relative_range       = 0.0;
  d->standard_deviation       = 0.0;
  d->coefficient_of_variation = 0.0;

  return d;
}

void ydata_analyze(ydata_t *d)
{
  if (!d)
    ynotif_error("ydata_analyze", YNOTIF_EXIT, "data pointer is NULL\n");
  
  _sort(d->datapoints, d->n_datapoints);

  d->min                      = d->datapoints[0];
  d->max                      = d->datapoints[d->n_datapoints - 1];
  d->mean                     = _mean(d->datapoints, d->n_datapoints);
  d->median                   = d->datapoints[(d->n_datapoints + 1) >> 1];
  d->variance                 = _variance(d->datapoints, d->n_datapoints, d->mean);
  d->min_median_diff          = fabs(d->min - d->median);
  d->min_relative_range       = (d->max - d->min) / d->min;
  d->standard_deviation       = sqrt(d->variance);
  d->coefficient_of_variation = (d->standard_deviation / d->mean) * 100.0;
}

void ydata_print(ydata_t *d)
{
  if (!d)
    ynotif_error("ydata_print", YNOTIF_EXIT, "data pointer is NULL\n");
  
  printf("%s; %lf; %lf; %lf; %lf; %lf; %lf; %lf; %lf; %lf\n",
	 (d->title) ? d->title : "",
	 d->min,
	 d->max,
	 d->mean,
	 d->median,
	 d->variance,
	 d->min_median_diff,
	 d->min_relative_range,
	 d->standard_deviation,
	 d->coefficient_of_variation);
}

void ydata_dump(ydata_t *d, const ascii *fname)
{
  if (!d)
    ynotif_error("ydata_dump", YNOTIF_EXIT, "data pointer is NULL\n");

  if (!fname)
    ynotif_error("ydata_dump", YNOTIF_EXIT, "file name pointer is NULL\n");

  FILE *fp = fopen(fname, "wb");

  if (!fp)
    ynotif_error("ydata_dump", YNOTIF_EXIT, "cannot create file '%s'\n", fname);
  
  fprintf(fp, "%s; %e; %e; %e; %e; %e; %e; %e; %e; %e\n",
	  (d->title) ? d->title : "",
	  d->min,
	  d->max,
	  d->mean,
	  d->median,
	  d->variance,
	  d->min_median_diff,
	  d->min_relative_range,
	  d->standard_deviation,
	  d->coefficient_of_variation);
}

f64 ydata_pearson_correlation(ydata_t *d1, ydata_t *d2)
{
  if (!d1)
    ynotif_error("ydata_pearson_correlation", YNOTIF_EXIT, "first data pointer is NULL\n");

  if (!d2)
    ynotif_error("ydata_pearson_correlation", YNOTIF_EXIT, "second data pointer is NULL\n");

  if (d1->n_datapoints != d2->n_datapoints)
    ynotif_error("ydata_pearson_correlation", YNOTIF_EXIT, "number of data points does not match\n");
  
  return _pearson(d1->datapoints, d2->datapoints, d1->n_datapoints);
}

void ydata_destroy(ydata_t **d)
{
  if (d)
    if (*d)
      {
	if ((*d)->datapoints)
	  free((*d)->datapoints);

	free(*d);
      }
    else
      ynotif_error("ydata_destroy", YNOTIF_EXIT, "data pointer is NULL\n");
  else
    ynotif_error("ydata_destroy", YNOTIF_EXIT, "pointer to data pointer is NULL\n");
}
