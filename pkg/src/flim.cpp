#include <RcppGSL.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>
#include <vector>

class Flim {
  // Actual parameters.
  RcppGSL::matrix<float> lambda_;
  RcppGSL::vector<float> kappa_;

  // Useful constants.
  RcppGSL::vector<float> ones_;

  // Outputs.
  RcppGSL::matrix<float> estimates_;

  // Temporaries.
  RcppGSL::vector<float> q_lambda_;
  RcppGSL::vector<float> singleton_expectation_;

  // Parameters.
  double beta1_;
  double beta2_;

  // Empirical data.
  RcppGSL::matrix<float> empirical_pair_;
  RcppGSL::vector<float> empirical_singleton_;

public:
  Flim(unsigned int N,
       double beta1,
       double beta2) : lambda_(N, N),
                       ones_(N),
                       kappa_(N),
                       estimates_(N, N),
                       q_lambda_(N),
                       singleton_expectation_(N),
                       beta1_(beta1),
                       beta2_(beta2),
                       empirical_pair_(N, N),
                       empirical_singleton_(N) {
    gsl_matrix_float_set_zero(lambda_);
    gsl_vector_float_set_zero(ones_);
  }

  ~Flim() {
    lambda_.free();
    ones_.free();
    kappa_.free();
    estimates_.free();
    q_lambda_.free();
    singleton_expectation_.free();
    empirical_pair_.free();
    empirical_singleton_.free();
  }
  
  // We need to compute something like:
  // sum_{k != i, j} \lambda_{i,k} x_i q_k(x_k = 1)
  // = x_i (\sum \lambda_{i, k} q_k(x_k = 1) - 
  //        \lambda_{i, j} q_j(x_j = 0))
  // (note that lambda_{i, i} = 0)
  unsigned int estimateExpectations() {
    // estimates_{x,y} = lambda_{x,y}
    gsl_matrix_float_memcpy(estimates_, lambda_);
    // estimates_{x,y} = lambda_{x,y} + kappa_x
    gsl_blas_sger(1.0, ones_, kappa_, estimates_);
    // estimates_{x,y} = lambda_{x,y} + kappa_x + kappa_y
    gsl_blas_sger(1.0, kappa_, ones_, estimates_);

    gsl_vector_float_set_zero(q_lambda_);
    gsl_blas_sgemv(CblasNoTrans,
                   1.0,
                   lambda_,
                   singleton_expectation_,
                   1.0,
                   q_lambda_);

    // estimates_{x,y} = lambda_{x,y} + kappa_x + kappa_y
    //                 + E_q [ sum lambda_{x, z} + lambda_{z, y} ]
    gsl_blas_sger(1.0, ones_, q_lambda_, estimates_);
    gsl_blas_sger(1.0, q_lambda_, ones_, estimates_);
  }

  void initializeKappa(int num_documents) {
    gsl_vector_float_memcpy(kappa_, empirical_singleton_);
    gsl_vector_float_scale(kappa_, num_documents);
    gsl_vector_float_add_constant(kappa_, 1.0);
    gsl_vector_float_scale(kappa_, 1.0 / (2.0 + num_documents));
    for (int ii = 0; ii < kappa_.size(); ++ii) {
      singleton_expectation_[ii] = kappa_[ii];
      kappa_[ii] = logit(kappa_[ii]);
    }
  }

  void loadCorpus(const std::vector<double>& singleton, 
                  const std::vector<double>& pair_x,
                  const std::vector<double>& pair_y,
                  const std::vector<double>& pair_count,
                  int num_documents) {
    for (int ii = 0; ii < singleton.size(); ++ii) {
      empirical_singleton_[ii] = singleton[ii] / num_documents;
    }
    for (int ii = 0; ii < pair_x.size(); ++ii) {
      int xx = pair_x[ii] - 1;
      int yy = pair_y[ii] - 1;
      empirical_pair_(xx, yy) = pair_count[ii] / num_documents;
    }
    initializeKappa(num_documents);
  }


  float sigmoid(float x) {
    return 1.0 / (1.0 + exp(-x));
  }

  float logit(float x) {
    return log(x / (1 - x));
  }

  // Normalizes precomputed statistics to give an expected value.
  float getComputedExpectation(unsigned int x, unsigned int y) {
    double overcount_x = lambda_(x,y) * singleton_expectation_[y];
    double overcount_y = lambda_(x,y) * singleton_expectation_[x];
    
    double p10 = q_lambda_[x] + kappa_[x] - overcount_x;
    double p01 = q_lambda_[y] + kappa_[y] - overcount_y;

    double p11 = estimates_(x,y) - overcount_x - overcount_y;

    return exp(p11) / (exp(p11) + exp(p10) + exp(p01) + 1);
  }

  void optimizeAll() {
    for (int x = 0; x < lambda_.nrow(); ++x) {
      for (int y = 0; y < x; ++y) {
        optimizeLambda(y, x);
      }
    }
  }

  void optimizeLambda(unsigned int x, unsigned int y) {
    double A = getComputedExpectation(x, y) * exp(-lambda_(x, y));
    double B = 2 * beta2_;
    double C = empirical_pair_(x, y) - beta1_;
    
    double delta1 = 0;
    double delta2 = 0;
    if (beta2_ == 0) {
      delta1 = log(C / A);
    } else {
      gsl_sf_result result;
      gsl_sf_lambert_W0_e(A / B * exp(C / B), &result);
      delta1 = C/B - result.val;
    }

    C += 2 * beta1_;
    if (beta2_ == 0) {
      delta2 = log(C / A);
    } else {
      gsl_sf_result result;
      gsl_sf_lambert_W0_e(A / B * exp(C / B), &result);
      delta2 = C/B - result.val;
    }

    double new_lambda = lambda_(x,y);
    if (new_lambda + delta1 >= 0) {
      new_lambda += delta1;
    } else if (new_lambda + delta2 <= 0) {
      new_lambda += delta2;
    } else {
      new_lambda = 0;
    }

    lambda_(x,y) = new_lambda;
    lambda_(y,x) = new_lambda;
  }

  RcppGSL::matrix<float> getLambda() {
    return lambda_;
  }
};

RCPP_MODULE(Rflim) {
	using namespace Rcpp;

	class_<Flim>("Flim")		
    .constructor<int,double,double>()
		.method("loadCorpus", &Flim::loadCorpus)
		.method("optimizeAll", &Flim::optimizeAll)
		.method("estimateExpectations", &Flim::estimateExpectations)
    .method("getLambda", &Flim::getLambda);
}
