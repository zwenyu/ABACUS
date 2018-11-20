#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List CptCinit(mat Y, int K, int numiter, double a_psi, double b_psi, List init = R_NilValue, bool useInit = false) {
  // Y is the p x N data matrix
  // K is the reduced dimension
  // d is the degree of differencing (up to 2)
  // burn_in, numiter for gibbs sampling
  // a_psi, b_psi are parameters for priors of psi 
  
  int p = Y.n_rows;
  int N = Y.n_cols;
  
  // make helper variables
  double sum1;
  double sum2;
  mat L(K,K);
  L.zeros();
  vec Me(K);
  Me.zeros();
  mat C_const(p,N);
  C_const.zeros();
  mat C_h(p,N);
  C_h.zeros();
  mat B_const(K,K);
  B_const.zeros();
  mat B_hinv(K,K);
  B_hinv.zeros();
  mat theta(p,N);
  theta.zeros();
  rowvec Ei(N);
  Ei.zeros();
  
  // Make and Initialize variables
  
  
  vec tau0(numiter);
  vec tau1(numiter);
  vec xi0(numiter);
  vec xi1(numiter);
  mat lam0(K, numiter);
  mat lam1(K, numiter);
  mat v0(K, numiter);
  mat v1(K, numiter);
  mat phi0(N,numiter);
  mat phi1(N,numiter);
  mat w0(N,numiter);
  mat w1(N,numiter);
  cube gamma0(N,K,numiter);
  cube gamma1(N,K,numiter);
  cube alpha0(N,K,numiter);
  cube alpha1(N,K,numiter);
  cube M(p,K,numiter);
  cube S(N,K,numiter);
  cube V0(N,K,numiter);
  cube V1(N,K,numiter);
  mat psi(p,numiter);
  
  tau0[0] = 1;
  tau1[0] = 1;
  xi0[0] = 1;
  xi1[0] = 1;
  lam0.col(0).ones();
  lam1.col(0).ones();
  v0.col(0).ones();
  v1.col(0).ones();
  phi0.col(0).ones();
  phi1.col(0).ones();
  w0.col(0).ones();
  w1.col(0).ones();
  gamma0.slice(0).ones();
  gamma1.slice(0).ones();
  alpha0.slice(0).ones();
  alpha1.slice(0).ones();
  M.slice(0).randn();
  S.slice(0).zeros();
  V0.slice(0).zeros();
  V1.slice(0).zeros();
  psi.col(0).ones();
  
  
  // Make difference operator
  
  mat D0(N,N);
  D0.eye();
  mat D1(N,N);
  D1.eye();
  
  vec subdiag(N-1);
  subdiag.fill(-1);
  D1.diag(-1) = subdiag;
  // D1 = D1 + diagmat(subdiag,-1);

  mat D1inv = D1.i();
  
  // Initialize Variables
  
  if(useInit){
    tau0[0] = as<double>(init["tau"]);
    xi0[0] = as<double>(init["xi"]);
    tau1[0] = as<double>(init["tau"]);
    xi1[0] = as<double>(init["xi"]);
    
    lam0.col(0) = as<vec>(init["lam"]);
    v0.col(0) = as<vec>(init["v"]);
    lam1.col(0) = as<vec>(init["lam"]);
    v1.col(0) = as<vec>(init["v"]);
    
    phi0.col(0) = as<vec>(init["phi0"]);
    w0.col(0) = as<vec>(init["w0"]);
    phi1.col(0) = as<vec>(init["phi1"]);
    w1.col(0) = as<vec>(init["w1"]);
    
    gamma0.slice(0) = as<mat>(init["gamma0"]);
    alpha0.slice(0) = as<mat>(init["alpha0"]);
    gamma1.slice(0) = as<mat>(init["gamma1"]);
    alpha1.slice(0) = as<mat>(init["alpha1"]);
    
    M.slice(0) = as<mat>(init["M"]);
    V0.slice(0) = as<mat>(init["V0"]);
    V1.slice(0) = as<mat>(init["V1"]);
    S.slice(0) = as<mat>(init["S"]);
    
    psi.col(0) = as<vec>(init["psi"]);
  }
  
  
  
  // Iterate
  
  for(int g = 0; g < numiter - 1; g++){
    
    // Fit xi, v, w
    
    xi0[g+1] = 1.0/R::rgamma(1, 1.0/(1.0+1.0/tau0[g]));
    xi1[g+1] = 1.0/R::rgamma(1, 1.0/(1.0+1.0/tau1[g]));
    
    for(int k = 0; k < K; k++){
      v0(k,g+1) = 1.0/R::rgamma(1, 1.0/(1.0+1.0/lam0(k,g)));
      v1(k,g+1) = 1.0/R::rgamma(1, 1.0/(1.0+1.0/lam1(k,g)));
    }
    for(int i = 0; i < N; i++){
      w0(i,g+1) = 1.0/R::rgamma(1, 1.0/(1.0+1.0/phi0(i,g)));
      w1(i,g+1) = 1.0/R::rgamma(1, 1.0/(1.0+1.0/phi1(i,g)));
    }
    
    // Fit tau
    
      // Fit tau0
      
    sum1 = 0;
    for(int i = 0; i < p; i++){
      for(int j = 0; j < K; j++){
        sum1 = sum1 + pow(M(i,j,g),2)/(2.0*psi(i,g) * lam0(j,g) * lam1(j,g) * tau1[g]);
      }
    } 
    
    sum2 = 0;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < K; j++){
        sum2 = sum2 + pow(V0(i,j,g),2)/(2.0*phi0(i,g)*lam0(j,g)*gamma0(i,j,g));
      }
    }
    
    tau0[g+1] = 1.0/R::rgamma((1+K*(p+N))/2.0, 1.0/(1.0/xi0[g+1] + sum1 + sum2));
    
      // Fit tau1
    sum1 = 0;
    for(int i = 0; i < p; i++){
      for(int j = 0; j < K; j++){
        sum1 = sum1 + pow(M(i,j,g),2)/(2.0*psi(i,g) * lam0(j,g) * lam1(j,g) * tau0[g+1]);
      }
    } 
    
    sum2 = 0;
    for(int i = 0; i < N; i++){
      for(int j = 0; j < K; j++){
        sum2 = sum2 + pow(V1(i,j,g),2)/(2.0*phi1(i,g)*lam1(j,g)*gamma1(i,j,g));
      }
    }
    
    tau1[g+1] = 1.0/R::rgamma((1+K*(p+N))/2.0, 1.0/(1.0/xi1[g+1] + sum1 + sum2));
    
    // Fit lam
    
      // Fit lam0
    for(int k = 0; k < K; k++){
      sum1 = 0;
      for(int i = 0; i < p; i++){
        sum1 = sum1 + pow(M(i,k,g),2)/(2*tau0[g+1]*tau1[g+1]*psi(i,g)*lam1(k,g));
      }
      sum2 = 0;
      for(int i = 0; i < N; i++){
        sum2 = sum2 + pow(V0(i,k,g),2)/(2*tau0[g+1]*phi0(i,g)*gamma0(i,k,g));
      }
      lam0(k,g+1) = 1.0/R::rgamma((1+p+N)/2, 1.0/(1.0/v0(k,g+1) + sum1 + sum2));
    }
      // Fit lam1
    for(int k = 0; k < K; k++){
      sum1 = 0;
      for(int i = 0; i < p; i++){
        sum1 = sum1 + pow(M(i,k,g),2)/(2*tau0[g+1]*tau1[g+1]*psi(i,g)*lam0(k,g+1));
      }
      sum2 = 0;
      for(int i = 0; i < N; i++){
        sum2 = sum2 + pow(V1(i,k,g),2)/(2*tau1[g+1]*phi1(i,g)*gamma1(i,k,g));
      }
      lam1(k,g+1) = 1.0/R::rgamma((1+p+N)/2.0, 1.0/(1.0/v1(k,g+1) + sum1 + sum2));
    }
    
    // Fit phi
      // Fit phi0
    for(int i = 0; i < N; i++){
      sum2 = 0;
      for(int k = 0; k < K; k++){
        sum2 = sum2 + pow(V0(i,k,g),2)/(2*tau0[g+1]*lam0(k,g+1)*gamma0(i,k,g));
      }
      phi0(i,g+1) = 1.0/R::rgamma((1+K)/2.0, 1.0/(1.0/w0(i,g+1) + sum2));
    }
      // Fit phi1
    for(int i = 0; i < N; i++){
      sum2 = 0;
      for(int k = 0; k < K; k++){
        sum2 = sum2 + pow(V1(i,k,g),2)/(2*tau1[g+1]*lam1(k,g+1)*gamma1(i,k,g));
      }
      phi1(i,g+1) = 1.0/R::rgamma((1+K)/2.0, 1.0/(1.0/w1(i,g+1) + sum2));
    }
    
    
    
    // Fit alpha
    
    for(int i = 0; i < N; i++){
      for(int k = 0; k < K; k++){
        alpha0(i,k,g+1) = 1.0/R::rgamma(1.0, 1.0/(1.0 + 1.0/gamma0(i,k,g)));
        alpha1(i,k,g+1) = 1.0/R::rgamma(1.0, 1.0/(1.0 + 1.0/gamma1(i,k,g)));
      }
    }
    
    
    // Fit gamma
    
    for(int i = 0; i < N; i++){
      for(int k = 0; k < K; k++){
        sum1 = 1.0/alpha0(i,k,g+1) + pow(V0(i,k,g),2)/(2*phi0(i,g+1)*lam0(k,g+1)*tau0[g+1]); 
        gamma0(i,k,g+1) = 1.0/R::rgamma(1.0, 1.0/sum1);
        sum1 = 1.0/alpha1(i,k,g+1) + pow(V1(i,k,g),2)/(2*phi1(i,g+1)*lam1(k,g+1)*tau1[g+1]); 
        gamma1(i,k,g+1) = 1.0/R::rgamma(1.0, 1.0/sum1);
      }
    }
    
    
    // Fit M
    
    mat QL(K,K);
    QL.zeros();
    QL.diag() = lam0.col(g+1) % lam1.col(g+1) * tau0[g+1] * tau1[g+1];
    L = (S.slice(g).t() * S.slice(g) + QL.i()).i(); 
    // L = (S.slice(g).t() * S.slice(g) + diagmat(lam0.col(g+1) % lam1.col(g+1) * tau0[g+1] * tau1[g+1]).i()).i(); 
    for(int i = 0; i < p; i++){
      Me = L*(S.slice(g).t() * Y.row(i).t());
      M.slice(g+1).row(i) = (Me + sqrt(psi(i,g)) * (chol(L).t() * randn(K))).t();
    }
    
    
    // Fit V
    
    C_const = Y - M.slice(g+1)*V0.slice(g).t() - M.slice(g+1)*V1.slice(g).t()*D1inv.t();
    mat QB(p,p);
    QB.zeros();
    QB.diag() = psi.col(g);
    B_const = M.slice(g+1).t() * QB.i() * M.slice(g+1);
    // B_const = M.slice(g+1).t() * diagmat(psi.col(g)).i() * M.slice(g+1);
    mat QB2(K,K);
    QB2.zeros();
    
    for(int i = 0; i < N; i++){
      // V1
      C_h = C_const + M.slice(g+1) * V1.slice(g).row(i).t() * D1inv.col(i).t();
      sum1 = 0;
      for(int j = 0; j < N; j++){
        sum1 = sum1 + pow(D1inv(j,i),2);
      }
      QB2.diag() = phi1(i,g+1) * tau1[g+1] * lam1.col(g+1) % gamma1.slice(g+1).row(i).t();
      B_hinv = (B_const * sum1 + (QB2.i())).i(); 
      V1.slice(g+1).row(i) = (B_hinv * (M.slice(g+1).t() * (QB.i()*C_h*D1inv.col(i))) + chol(B_hinv).t() * randn(K)).t();
      // B_hinv = (B_const * sum1 + (diagmat(phi1(i,g+1) * tau1[g+1] * lam1.col(g+1) % gamma1.slice(g+1).row(i).t()).i())).i(); 
      // V1.slice(g+1).row(i) = (B_hinv * (M.slice(g+1).t() * (diagmat(psi.col(g)).i()*C_h*D1inv.col(i))) + chol(B_hinv).t() * randn(K)).t();
     
      // Update C_const
      C_const = C_h - M.slice(g+1) * V1.slice(g+1).row(i).t() * D1inv.col(i).t();
      
      // V0
      C_h = C_const + M.slice(g+1) * V0.slice(g).row(i).t() * D0.col(i).t();
      QB2.diag() = phi0(i,g+1) * tau0[g+1] * lam0.col(g+1) % gamma0.slice(g+1).row(i).t();
      B_hinv = (B_const + (QB2.i())).i(); 
      V0.slice(g+1).row(i) = (B_hinv * (M.slice(g+1).t() * (QB.i()*C_h*D0.col(i))) + chol(B_hinv).t() * randn(K)).t();
      // B_hinv = (B_const + (diagmat(phi0(i,g+1) * tau0[g+1] * lam0.col(g+1) % gamma0.slice(g+1).row(i).t()).i())).i(); 
      // V0.slice(g+1).row(i) = (B_hinv * (M.slice(g+1).t() * (diagmat(psi.col(g)).i()*C_h*D0.col(i))) + chol(B_hinv).t() * randn(K)).t();
      
      // Update C_const
      C_const = C_h - M.slice(g+1) * V0.slice(g+1).row(i).t() * D0.col(i).t();
    }
    
    // Fit S
    
    S.slice(g+1) = V0.slice(g+1) + D1inv*V1.slice(g+1);
    
    // Fit Psi
    
    theta = M.slice(g+1) * S.slice(g+1).t();
    
    for(int i = 0; i < p; i++){
      Ei = Y.row(i) - theta.row(i);
      psi(i,g+1) = 1.0/R::rgamma(a_psi + N/2, 1.0/(b_psi + .5 * dot(Ei,Ei)));
    }
    
    
  }
  
  List outList;
  outList["tau0"] = tau0;
  outList["tau1"] = tau1;
  outList["xi0"] = xi0;
  outList["xi1"] = xi1;
  outList["lam0"] = lam0;
  outList["lam1"] = lam1;
  outList["v0"] = v0;
  outList["v1"] = v1;
  outList["phi0"] = phi0;
  outList["phi1"] = phi1;
  outList["w0"] = w0;
  outList["w1"] = w1;
  outList["gamma0"] = gamma0;
  outList["gamma1"] = gamma1;
  outList["alpha0"] = alpha0;
  outList["alpha1"] = alpha1;
  outList["M"] = M;
  outList["V0"] = V0;
  outList["V1"] = V1;
  outList["S"] = S;
  outList["psi"] = psi;
  
  return outList;
  
  
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


/*** R
*/
