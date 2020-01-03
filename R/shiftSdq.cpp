#include <iostream>
#include <math.h>
#include <bits/stdc++.h>
#include <rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix shiftSdq(const NumericMatrix& m1){
	NumericMatrix out(m1.nrow(),m1.ncol());
  
	for(double i=0; i <= (m1.nrow() - 1); i++) {
	   	for(double j=0; j <= (m1.ncol() - 1); j++) {
	   		if (std::isnan(m1(i,j))) {
	   			out(i,j) = std::nan("");
        } else {
		        out(i,j) = pow((m1(i,j) - m1(i-1,j)), 2) + pow((m1(i,j) - m1(i,j-1)), 2);
		    }
		}
	}
	return out;
}