/*********************************************************************************
 *
 * burns.c - Defines a bunch of burns that you can apply to initial conditions
 * Copyright(c) 2014 Kevin Croker
 * GPL v3
 *
 **********************************************************************************/

#include <math.h>
#include "tkgadget2.h"
#define PI32 5.56832799 

float gaussianBurn(float radius, float *params) {
  
  // gaussian 
  // <sigma> <amplitude>
  // return 1.0/PI32 * exp(-0.5 * gsl_pow_2(radius / atof(argv[0]))) * atof(argv[1]);
  return 0.0;
}

float tophatBurn(float radius, float *params) {
  
  char where;
  float chop;

  chop = params[0];
  where = params[1] > 0;
  
  // Chop inside
  if(where) {
    if(radius < chop)
      return 0.0;
    else
      return 1.0;
  }
  else {

    if(radius > chop)
      return 0.0;
    else
      return 1.0;
  }
}

float thinBurn(float radius, float *params) {

  // Thins a region outside params[0] by a factor params[1]
  if(radius >= params[0])
    return params[1];
  else
    return 1.0;
}
