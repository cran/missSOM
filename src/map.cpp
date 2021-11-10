/* map.c: calculate distances of objects to units, now using the
   appropriate distance functions. We also need to take into account
   weights. We do not return all distances but only the winners and
   the associated overall distances.

   This version: July 27, 2016
   Author: Ron Wehrens and Johannes Kruisselbrink
*/

#include <Rcpp.h>

#include "missSOM.h"
#include "distance-functions.h"

Rcpp::List RcppMap(
    Rcpp::NumericMatrix data,   /* objects to be mapped */
    Rcpp::NumericMatrix codes,
    Rcpp::ExpressionVector distanceFunctions)
{
  int
    numObjects = data.ncol(),     /* number of objects */
    numCodes = codes.ncol(),      /* number of units in the map */
    totalVars = data.nrow(),      /* total number of variables sum(numVars) */
    i, nearest;

  double distance;

  Rcpp::IntegerVector winners(numObjects);
  Rcpp::NumericVector unitDistances(numObjects);

  double *pCodes = REAL(codes);


  /* Get the distance function pointers. */
  std::vector<DistanceFunctionPtr> distanceFunctionPtrs =
    GetDistanceFunctions(distanceFunctions);


  /* Loop over all data objects */
  for (i = 0; i < numObjects; i++) {
    /* Find best matching unit index and distance */
    FindBestMatchingUnit(
      &data[i * totalVars],
      pCodes,
      numCodes,
      totalVars,
      distanceFunctionPtrs,
      nearest,
      distance);

    winners[i] = nearest;
    unitDistances[i] = distance;
  }

  return Rcpp::List::create(
    Rcpp::Named("winners") = winners,
    Rcpp::Named("unitdistances") = unitDistances);
}
