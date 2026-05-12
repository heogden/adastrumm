// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math/mix.hpp> // need this to be able to do hessian
#include <stan/math.hpp>  // pulls in everything from rev/ and prim/
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> find_u(const Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha) {
  T alpha_norm = stan::math::sqrt(stan::math::dot_self(alpha));
  Eigen::Matrix<T, Eigen::Dynamic, 1> u = alpha;
  u(0) -= alpha_norm;
  return u;
}

template <typename T>
T find_gamma(const Eigen::Matrix<T, Eigen::Dynamic, 1>& u) {
  T t = stan::math::dot_self(u);
  T gamma = 2 / t;
  return gamma;
}

template <typename T> struct HouseholderReduced {
  // members
  const size_t n_rows; // the matrix version has n_rows rows and n_rows -1 columns
  const Eigen::Matrix<T, Eigen::Dynamic, 1> u;
  const T gamma;
  
  // constructor
  HouseholderReduced(const Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha): n_rows(alpha.size()),
								  u(find_u(alpha)),
								  gamma(find_gamma(u))
  {
  }

  // multiply by a vector x
  Eigen::Matrix<T, Eigen::Dynamic, 1> operator*(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> x_ext(x.size() + 1);
    x_ext << 0, x;
    T a = u.dot(x_ext);
    return x_ext - a * gamma * u;
  }
  
};


template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> find_beta(const Eigen::Matrix<T, Eigen::Dynamic, 1>& alpha,
					      size_t K, size_t n_B) {
  std::vector<HouseholderReduced<T> > Hstar_vec;
  Eigen::Matrix<T, Eigen::Dynamic, 1> beta(K * n_B);
  
  for(size_t k = 0; k < K; ++k) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> alpha_k = alpha.segment(k * n_B - k * (k-1) / 2, n_B - k);
    
    Eigen::Matrix<T, Eigen::Dynamic, 1> beta_k = alpha_k;
    for(auto it = Hstar_vec.rbegin(); it != Hstar_vec.rend(); ++it) {
      beta_k = (*it) * beta_k;
    }

    beta.segment(n_B * k, n_B) = beta_k;

    if (k < K - 1) {
      HouseholderReduced<T> Hstar_k(alpha_k);
      Hstar_vec.push_back(Hstar_k);
    }
  
  }
  return beta;
}



struct cluster {
  // members
  size_t n_c;
  Eigen::MatrixXd X_c;
  Eigen::VectorXd y_c;

  // constuctor
  cluster() {
    n_c = 0;
  }

  void append(const Eigen::RowVectorXd& X_i, double y_i, size_t n_B) {  
    n_c += 1;

    X_c.conservativeResize(n_c, n_B);
    X_c.row(n_c - 1) = X_i;
    y_c.conservativeResize(n_c, 1);
    y_c(n_c - 1) = y_i;
  }
  
};


struct loglikp_func {
  // members
  const size_t n_B;
  std::vector<cluster> clusters;
  double sp;
  Eigen::MatrixXd S;
  const size_t K;

  // constructor
  loglikp_func(Eigen::MatrixXd& X, Eigen::VectorXd& y, std::vector<int>& c, double sp_, Eigen::MatrixXd& S_, size_t K_): n_B(X.cols()), sp(sp_), S(S_), K(K_) {
    auto d = max_element(c.begin(), c.end());
    clusters.resize(*d + 1);
    for(size_t i = 0; i < c.size(); ++i) {
      int cluster_id = c[i];
      clusters[cluster_id].append(X.row(i), y(i), n_B);
    }
  }

  // function definition
  template <typename T>
  T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta) const {
    T log_sigma = theta.tail(1)[0];
    T sigma = stan::math::exp(log_sigma);
    T tau= 1 / (sigma * sigma);

    Eigen::Matrix<T, Eigen::Dynamic, 1> beta_0 = theta.head(n_B);

    // The correct number of alpha parameters is K * n_B - K * (K - 1) / 2
    Eigen::Matrix<T, Eigen::Dynamic, 1> alpha = theta.segment(n_B, K * n_B - K * (K - 1) / 2);
    Eigen::Matrix<T, Eigen::Dynamic, 1> beta = find_beta(alpha, K, n_B);
    
    std::vector<T> l_contribs;
    
    for(auto it = clusters.begin(); it < clusters.end(); it++) {
      cluster clust = *it;
      Eigen::Matrix<T, Eigen::Dynamic, 1> f0_c = clust.X_c * beta_0;
      Eigen::Matrix<T, Eigen::Dynamic, 1> z = clust.y_c - f0_c;

      T ldet_Sigma = 2.0 * clust.n_c * log_sigma;
      T Q = tau * z.dot(z);
    
      std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> D;
    
      for(size_t k = 0; k < K; ++k) {
	Eigen::Matrix<T, Eigen::Dynamic, 1> beta_k = beta.segment(n_B * k, n_B);
	Eigen::Matrix<T, Eigen::Dynamic, 1> f_k = clust.X_c * beta_k;
	Eigen::Matrix<T, Eigen::Dynamic, 1> b_k = tau * f_k;
	for(int j = 0; j < k; ++j) {
	  Eigen::Matrix<T, Eigen::Dynamic, 1> d_j = D[j];
	  T d_j_f_k = d_j.dot(f_k); 
	  b_k = b_k - d_j_f_k * d_j;
	}
	T a_k = 1 + b_k.dot(f_k);

        // Numerical stability check for a_k
        if (a_k <= 0) {
          // Return -INFINITY for log-likelihood if a_k is not positive definite
          return stan::math::negative_infinity();
        }

	Eigen::Matrix<T, Eigen::Dynamic, 1> d_k = b_k / stan::math::sqrt(a_k);
	D.push_back(d_k);

	ldet_Sigma += log(a_k);
	T d_k_z = d_k.dot(z);
	Q -= d_k_z * d_k_z;
      }
    
      l_contribs.push_back(- (clust.n_c * stan::math::log(2 * M_PI) + ldet_Sigma + Q) / 2);
    }
    T l = stan::math::sum(l_contribs);
       
    std::vector<T> w_contribs;
    T w_0 =  stan::math::quad_form(S, beta_0);
    w_contribs.push_back(w_0);
    
    for(size_t k = 0; k < K; ++k) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> beta_k = beta.segment(n_B * k, n_B);
      
      T w_k = stan::math::quad_form(S, beta_k);
      w_contribs.push_back(w_k);
    }
    T w = stan::math::sum(w_contribs);
    T pen = 0.5 * sp * tau * w;
      
    return l - pen;
  }
};
  
// [[Rcpp::export]]
double loglikelihood_pen(Eigen::VectorXd theta,
			 Eigen::MatrixXd X,
			 Eigen::VectorXd y,
			 std::vector<int> c,
			 double sp,
			 Eigen::MatrixXd S,
			 size_t K) {

  loglikp_func lf(X, y, c, sp, S, K);
  return lf(theta);
}


// [[Rcpp::export]]
NumericVector loglikelihood_pen_grad(Eigen::VectorXd theta,
				     Eigen::MatrixXd X,
				     Eigen::VectorXd y,
				     std::vector<int> c,
				     double sp,
				     Eigen::MatrixXd S,
				     size_t K) {
  // declarations
  double l;
  Eigen::VectorXd grad_l;

  loglikp_func lf(X, y, c, sp, S, K);

  stan::math::gradient(lf, theta, l, grad_l);

  
  // reformat returned result
  NumericVector grad_l1 = wrap(grad_l);
  
  return grad_l1;
}


// [[Rcpp::export]]
NumericVector loglikelihood_pen_hess(Eigen::VectorXd theta,
				     Eigen::MatrixXd X,
				     Eigen::VectorXd y,
				     std::vector<int> c,
				     double sp,
				     Eigen::MatrixXd S,
				     size_t K) {
  // declarations
  double l;
  Eigen::VectorXd grad_l;
  Eigen::MatrixXd hess_l;

  loglikp_func lf(X, y, c, sp, S, K);

  stan::math::hessian(lf, theta, l, grad_l, hess_l);

  
  // reformat returned result
  NumericMatrix hess_l1 = wrap(hess_l);

  return hess_l1;
}
