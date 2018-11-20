#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List CptC(mat Y, int K, int d, int numiter, double a_psi, double b_psi) {
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
	
	vec tau(numiter);
	vec xi(numiter);
	mat lam(K, numiter);
	mat v(K, numiter);
	mat phi(N,numiter);
	mat w(N,numiter);
	cube gamma(N,K,numiter);
	cube alpha(N,K,numiter);
	cube M(p,K,numiter);
	cube S(N,K,numiter);
	cube V(N,K,numiter);
	mat psi(p,numiter);
	
	tau[0] = 1;
	xi[0] = 1;
	lam.col(0).ones();
	v.col(0).ones();
	phi.col(0).ones();
	w.col(0).ones();
	gamma.slice(0).ones();
	alpha.slice(0).ones();
	M.slice(0).randn();
	S.slice(0).zeros();
	V.slice(0).zeros();
	psi.col(0).ones();

	// Make difference operator
	
	mat D(N,N);
	D.eye();
	
	if(d == 1){
		vec subdiag(N-1);
		subdiag.fill(-1);
		D.diag(-1) -= 1;
	}else if(d==2){
		vec subdiag(N-1);
		vec subsubdiag(N-2);
		subdiag.fill(-2);
		subsubdiag.fill(1);
		D.diag(-1) -= 2;
		D.diag(-2) += 1;
	}
	mat Dinv = D.i();
	
	// Iterate
	
	for(int g = 0; g < numiter - 1; g++){
		
		// Fit xi, v, w
		
		xi[g+1] = 1.0/R::rgamma(1, 1.0/(1.0+1.0/tau[g]));
		for(int k = 0; k < K; k++){
			v(k,g+1) = 1.0/R::rgamma(1, 1.0/(1.0+1.0/lam(k,g)));
		}
		for(int i = 0; i < N; i++){
			w(i,g+1) = 1.0/R::rgamma(1, 1.0/(1.0+1.0/phi(i,g)));
		}

		// Fit tau
		
		sum1 = 0;
		for(int i = 0; i < p; i++){
			for(int j = 0; j < K; j++){
				sum1 = sum1 + pow(M(i,j,g),2)/(2.0*psi(i,g) * lam(j,g));
			}
		} 
		
		sum2 = 0;
		for(int i = 0; i < N; i++){
			for(int j = 0; j < K; j++){
				sum2 = sum2 + pow(V(i,j,g),2)/(2.0*phi(i,g)*lam(j,g)*gamma(i,j,g));
			}
		}
		
		tau[g+1] = 1.0/R::rgamma((1+K*(p+N))/2.0, 1.0/(1.0/xi[g+1] + sum1 + sum2));
		
		// Fit lam
		
		for(int k = 0; k < K; k++){
			sum1 = 0;
			for(int i = 0; i < p; i++){
				sum1 = sum1 + pow(M(i,k,g),2)/(2*tau[g+1]*psi(i,g));
			}
			sum2 = 0;
			for(int i = 0; i < N; i++){
				sum2 = sum2 + pow(V(i,k,g),2)/(2*tau[g+1]*phi(i,g)*gamma(i,k,g));
			}
			lam(k,g+1) = 1.0/R::rgamma((1+p+N)/2.0, 1.0/(1.0/v(k,g+1) + sum1 + sum2));
		}
		
		// Fit phi
		
		for(int i = 0; i < N; i++){
			sum2 = 0;
			for(int k = 0; k < K; k++){
				sum2 = sum2 + pow(V(i,k,g),2)/(2*tau[g+1]*lam(k,g+1)*gamma(i,k,g));
			}
			phi(i,g+1) = 1.0/R::rgamma((1+K)/2.0, 1.0/(1.0/w(i,g+1) + sum2));
		}
		
		// Fit alpha
		
		for(int i = 0; i < N; i++){
			for(int k = 0; k < K; k++){
				alpha(i,k,g+1) = 1.0/R::rgamma(1.0, 1.0/(1.0 + 1.0/gamma(i,k,g)));
			}
		}
		
		// Fit gamma
		
		for(int i = 0; i < N; i++){
			for(int k = 0; k < K; k++){
				sum1 = 1.0/alpha(i,k,g+1) + pow(V(i,k,g),2)/(2*phi(i,g+1)*lam(k,g+1)*tau[g+1]); 
				gamma(i,k,g+1) = 1.0/R::rgamma(1.0, 1.0/sum1);
			}
		}
		
		// Fit M
		
		mat QL(K,K);
		QL.zeros();
		QL.diag() = lam.col(g+1) * tau[g+1];
		L = (S.slice(g).t() * S.slice(g) + QL.i()).i(); 
		for(int i = 0; i < p; i++){
			Me = L*(S.slice(g).t() * Y.row(i).t());
			M.slice(g+1).row(i) = (Me + sqrt(psi(i,g)) * (chol(L).t() * randn(K))).t();
		}		
		
		// Fit V
		
		C_const = Y - M.slice(g+1)*V.slice(g).t()*Dinv.t();
		mat QB(p,p);
		QB.zeros();
		QB.diag() = psi.col(g);
		B_const = M.slice(g+1).t() * QB.i() * M.slice(g+1);
		mat QB2(K,K);
		QB2.zeros();
		
		for(int i = 0; i < N; i++){
			C_h = C_const + M.slice(g+1) * V.slice(g).row(i).t() * Dinv.col(i).t();
			sum1 = 0;
			for(int j = 0; j < N; j++){
				sum1 = sum1 + pow(Dinv(j,i),2);
			}
			QB2.diag() = phi(i,g+1) * tau[g+1] * lam.col(g+1) % gamma.slice(g+1).row(i).t();
			B_hinv = (B_const * sum1 + (QB2.i())).i(); 
			V.slice(g+1).row(i) = (B_hinv * (M.slice(g+1).t() * (QB.i()*C_h*Dinv.col(i))) + chol(B_hinv).t() * randn(K)).t();

			// Update C_const
			C_const = C_h - M.slice(g+1) * V.slice(g+1).row(i).t() * Dinv.col(i).t();
		}
		
		// Fit S
		
		S.slice(g+1) = Dinv*V.slice(g+1);
		
		// Fit Psi
		
		theta = M.slice(g+1) * S.slice(g+1).t();
		
		for(int i = 0; i < p; i++){
			Ei = Y.row(i) - theta.row(i);
			psi(i,g+1) = 1.0/R::rgamma(a_psi + N/2.0, 1.0/(b_psi + .5 * dot(Ei,Ei)));
		}
		
	}
	
	return Rcpp::List::create(Rcpp::Named("tau") = tau,
		Rcpp::Named("xi") = xi,
		Rcpp::Named("lam") = lam,
		Rcpp::Named("v") = v,
		Rcpp::Named("phi") = phi,
		Rcpp::Named("w") = w,
		Rcpp::Named("gamma") = gamma,
		Rcpp::Named("alpha") = alpha,
		Rcpp::Named("M") = M,
		Rcpp::Named("V") = V,
		Rcpp::Named("S") = S,
		Rcpp::Named("psi") = psi);
		
}


/*** R
*/
