#include <TMB.hpp>
using namespace density;
template<class Type>
Type objective_function<Type>::operator() () {
  //------
  // DATA
  //------
  DATA_MATRIX(y); // data
  DATA_MATRIX(ypresent); // observations present
  DATA_IVECTOR(first_obs); // indicator for first observations
  //------------
  // PARAMETERS 
  //------------
  PARAMETER(lnsde);
  PARAMETER(lnsdx);
  PARAMETER(logitrho);
  PARAMETER_VECTOR(x); // underlying level
  PARAMETER_VECTOR(Apar); // stock fixed effects
  //--------------------------
  // PRELIMINARY CALCULATIONS
  //--------------------------
  int n = y.cols();
  int m = y.rows();
  Type sde = exp(lnsde);
  Type sdx = exp(lnsdx);
  Type rho = -1.0 + 2.0 / (1.0 + exp(-logitrho)); // bounded [-1,1]
  // A parameters
  vector<Type> A(m);
  //A << Type(0), Apar;
  A << (-sum(Apar)), Apar; // note sum-to-zero constraint
  //-----------
  // PROCEDURE 
  //-----------
  Type nll = 0.0; // initialize negative log likelihood
  // random walk trend
  for(int i = 1; i < n; i++){
    nll -= dnorm(x(i), x(i-1), sdx, true);
  }
  // observations
  for(int j = 0; j < m; j++){
    // see Greene equation 14-63
    // first observation
    nll -= -0.5 * (log(2.0 * M_PI) + log(pow(sde, 2.0)) - log(1 - pow(rho, 2.0))) - 
      (1.0 - pow(rho, 2.0)) / (2.0 * pow(sde, 2.0)) * pow((y(j,first_obs(j)) - (x(first_obs(j)) + A(j))), 2.0);
    // other observations
    for(int i = (first_obs(j) + 1); i < n; i++){
      if((ypresent(j,i) > 0 & ypresent(j,i-1) > 0)){
	nll -= -0.5 * (log(2.0 * M_PI) + log(pow(sde, 2.0))) -
	  1.0 / (2.0 * pow(sde, 2.0)) * pow(((y(j,i) - rho * y(j,i-1)) - (x(i) - rho * x(i-1)) - (A(j) - rho * A(j))), 2.0);
      }
    }
  }
  ADREPORT(A);
  return nll;
}
