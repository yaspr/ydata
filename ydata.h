#pragma once

#include "types.h"

#define YDATA_MAX_TITLE      64
#define YDATA_MIN_DATAPOINTS 33

typedef struct ydata_s {

  u64 n_datapoints;
  
  const ascii *title;
  
  f64 *datapoints;
  
  f64 min;
  f64 max;
  f64 mean;
  f64 median;
  f64 variance;
  f64 min_median_diff;
  f64 min_relative_range;
  f64 standard_deviation;
  f64 coefficient_of_variation;
  
} ydata_t;

ydata_t *ydata_create();
void ydata_analyze(ydata_t *d);
void ydata_print(ydata_t *d);
void ydata_dump(ydata_t *d, const ascii *fname);
f64  ydata_pearson_correlation(ydata_t *d1, ydata_t *d2);
void ydata_destroy(ydata_t **d);
