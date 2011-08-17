#include <Rcpp.h>

#include <RcppGSL.h>
#include <<gsl/gsl_matrix.h>

class Flim {
  RcppGSL::matrix<float> lambda_(N, N);
  RcppGSL::vector<float> ones_(N);

public:
  unsigned int estimateExpectations() {
      
  }
};

RCPP_MODULE(Rflim) {
	using namespace Rcpp;

	class_<Flim>("Flim")		
		.method("estimateExpectations", &Flim::estimateExpectations);
}
