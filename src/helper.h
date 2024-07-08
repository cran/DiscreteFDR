// //[[Rcpp::depends(RcppClock)]]
#include <RcppArmadillo.h>
//#include <RcppClock.h>
using namespace Rcpp;

// fast step function evaluation
// step function is represented by a single numeric vector under the conditions
// a) f(x) = x and b) x is sorted
// this is much faster than passing and evaluating R step function objects
// struct stepfun_result{
//   NumericVector y;
//   int last_index;
// };
// 
// stepfun_result stepfun(const NumericVector &x, const NumericVector &sfun, const int &start = 0);

// NumericVector stepfun(const NumericVector &x, const NumericVector &sfun);

// int stepfun_idx(const double &x, const NumericVector &sfun, const int &start)

// shortcut function that eliminates all values of a SORTED vector that
// are < limit, except the largest value <= limit
// NumericVector short_eff(const NumericVector &x, const double limit);

inline void eval_pv(double &eval, double val, const NumericVector &vec, int len, int &pos){
  //if(val < 1){
    while(pos < len && vec[pos] <= 1 && vec[pos] <= val) pos++;
    if(pos) eval = vec[pos - 1];
    else eval = 0;
  //} else eval = 1;
}

inline void eval_pv_rev(double &eval, double val, const NumericVector &vec, int &pos){
  //if(val < 1){
    while(pos > 0 && (vec[pos] > val || vec[pos] > 1)) pos--;
    if(vec[pos] <= val) eval = vec[pos];
    else eval = 0;
  //} else eval = 1;
}

// function that binds two vectors, sorts it and eliminates duplications 
NumericVector sort_combine(const NumericVector &x, const NumericVector &y);

// sort columns of a matrix in descending order
// using an intermediate numeric vector is necessary, because "in-column"
// sorting does not always work as expected, especially for large columns
//void colsortdec(NumericMatrix &mat);

// sort order
IntegerVector order(const NumericVector &x, bool descending = false);
