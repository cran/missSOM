// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "missSOM_types.h"
#include <Rcpp.h>

using namespace Rcpp;

// CreateStdDistancePointers
Rcpp::ExpressionVector CreateStdDistancePointers(const Rcpp::IntegerVector& types);
RcppExport SEXP _missSOM_CreateStdDistancePointers(SEXP typesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type types(typesSEXP);
    rcpp_result_gen = Rcpp::wrap(CreateStdDistancePointers(types));
    return rcpp_result_gen;
END_RCPP
}
// CreateStdDistancePointer
Rcpp::XPtr<DistanceFunctionPtr> CreateStdDistancePointer(int type);
RcppExport SEXP _missSOM_CreateStdDistancePointer(SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(CreateStdDistancePointer(type));
    return rcpp_result_gen;
END_RCPP
}
// ObjectDistances
Rcpp::NumericVector ObjectDistances(Rcpp::NumericMatrix data, Rcpp::ExpressionVector distanceFunctions);
RcppExport SEXP _missSOM_ObjectDistances(SEXP dataSEXP, SEXP distanceFunctionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::ExpressionVector >::type distanceFunctions(distanceFunctionsSEXP);
    rcpp_result_gen = Rcpp::wrap(ObjectDistances(data, distanceFunctions));
    return rcpp_result_gen;
END_RCPP
}
// RcppImputeSOM
Rcpp::List RcppImputeSOM(Rcpp::NumericMatrix data, Rcpp::NumericMatrix missData, Rcpp::NumericMatrix codes, Rcpp::ExpressionVector distanceFunctions, Rcpp::NumericMatrix neighbourhoodDistances, int neighbourhoodFct, Rcpp::NumericVector alphas, Rcpp::NumericVector radii, int numEpochs, bool bool_impute, Rcpp::IntegerVector missingCol, Rcpp::IntegerVector missingRow);
RcppExport SEXP _missSOM_RcppImputeSOM(SEXP dataSEXP, SEXP missDataSEXP, SEXP codesSEXP, SEXP distanceFunctionsSEXP, SEXP neighbourhoodDistancesSEXP, SEXP neighbourhoodFctSEXP, SEXP alphasSEXP, SEXP radiiSEXP, SEXP numEpochsSEXP, SEXP bool_imputeSEXP, SEXP missingColSEXP, SEXP missingRowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type missData(missDataSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type codes(codesSEXP);
    Rcpp::traits::input_parameter< Rcpp::ExpressionVector >::type distanceFunctions(distanceFunctionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type neighbourhoodDistances(neighbourhoodDistancesSEXP);
    Rcpp::traits::input_parameter< int >::type neighbourhoodFct(neighbourhoodFctSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type radii(radiiSEXP);
    Rcpp::traits::input_parameter< int >::type numEpochs(numEpochsSEXP);
    Rcpp::traits::input_parameter< bool >::type bool_impute(bool_imputeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type missingCol(missingColSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type missingRow(missingRowSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppImputeSOM(data, missData, codes, distanceFunctions, neighbourhoodDistances, neighbourhoodFct, alphas, radii, numEpochs, bool_impute, missingCol, missingRow));
    return rcpp_result_gen;
END_RCPP
}
// RcppMap
Rcpp::List RcppMap(Rcpp::NumericMatrix data, /* objects to be mapped */     Rcpp::NumericMatrix codes, Rcpp::ExpressionVector distanceFunctions);
RcppExport SEXP _missSOM_RcppMap(SEXP dataSEXP, SEXP codesSEXP, SEXP distanceFunctionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< /* objects to be mapped */     Rcpp::NumericMatrix >::type codes(codesSEXP);
    Rcpp::traits::input_parameter< Rcpp::ExpressionVector >::type distanceFunctions(distanceFunctionsSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppMap(data, codes, distanceFunctions));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_missSOM_CreateStdDistancePointers", (DL_FUNC) &_missSOM_CreateStdDistancePointers, 1},
    {"_missSOM_CreateStdDistancePointer", (DL_FUNC) &_missSOM_CreateStdDistancePointer, 1},
    {"_missSOM_ObjectDistances", (DL_FUNC) &_missSOM_ObjectDistances, 2},
    {"_missSOM_RcppImputeSOM", (DL_FUNC) &_missSOM_RcppImputeSOM, 12},
    {"_missSOM_RcppMap", (DL_FUNC) &_missSOM_RcppMap, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_missSOM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}