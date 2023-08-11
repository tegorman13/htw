#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double sumC(NumericVector x)
{
  
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}

// [[Rcpp::export]]
NumericVector pointVecDistCpp(double stimValue, NumericVector vec) {
  NumericVector out;
  out = vec - stimValue;  
  return out;
}




// [[Rcpp::export]]
NumericVector inputActivationCpp(double c, NumericVector inDist) {
  double sumAct; NumericVector distSq;
  NumericVector inputActivations(inDist.size());
  distSq = (-c) *(inDist*inDist);
  inputActivations = exp(distSq);
  sumAct = sumC(inputActivations);
  inputActivations = inputActivations/sumAct;
  return inputActivations;
}


//lookupBandBounds <- NumericVector(band)




  
  
// Calculate activation of output nodes
// [[Rcpp::export]]
NumericVector aocalcCpp(NumericMatrix w, NumericVector ah) {
  int i,j, nrow = w.nrow(), ncol = w.ncol();
  NumericVector out(nrow);
  for(j=0;j < nrow; j++) {
    out(j) = 0.0;
    for(i=0;i < ncol; i++) {
      out(j) += ah(i) * w(j,i);
    }
  }
  return out;
}

// [[Rcpp::export]]
double almCpp(NumericVector outputActivations, NumericVector outputNode)
{
  double sumAct; double almPred; NumericVector outProbs;
  NumericVector nodeOut;
  sumAct = sumC(outputActivations);
  outProbs = outputActivations/sumAct;
  nodeOut = outputNode*outProbs;
  //print(nodeOut);
  almPred = sumC(nodeOut);
  return(almPred);
}


// 
// double inputToALMCpp (double stimValue, double c, NumericVector vec,NumericVector outputNode, NumericMatrix w) 
// {
//   
//   NumericVector inp; NumericVector ia; NumericVector oa;
//   double aOut;
//   
//   inp = pointVecDistCpp(stimValue,vec);
//   ia = inputActivationCpp(c, inp);
//   oa = aocalcCpp(w, ia);
//   aOut = almCpp(oa,outputNode);
//   
//   return(aOut);
//   
// }
//   
  
  
  
  
  
  
  // [[Rcpp::export]]
  NumericMatrix deltawcalcCpp(double lw, NumericVector pe, NumericVector ah) {
    int i,j, nout = pe.size(), nhid = ah.size();
    NumericMatrix delt(nout,nhid);
    for(i=0; i < nout; i++) {
      for(j=0; j < nhid; j++) {
        delt(i,j) = lw * pe[i] * ah[j];
      }
    }
    return delt;
  }

  
  
  
  // Update associative link strength
  // [[Rcpp::export]]
  NumericMatrix wupdateCpp(NumericMatrix w, NumericMatrix delt) {
    int i, j, nhid = w.ncol(), nout = w.nrow();
    NumericMatrix nw(nout,nhid);
    for(i=0; i < nout; i++) {
      for(j=0; j < nhid; j++) {
        nw(i,j) = w(i,j) + delt(i,j);
      }
    }
    return nw;
  }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# sv = 10;
# iv <- seq(1,10)
# wm <- matrix(2,5,length(iv))
# outNodes = seq(1,10)
# d=inputDist(sv,iv)
# ia=inputActivation(.1,d)
# oa<- aocalc(wm,ia)
# oa/sumC(oa)
*/
