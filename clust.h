#include <gsl/gsl_rng.h>
#include "armadillo"
#include "imf_pflamm.h"

using namespace arma;

class clust
{
	public:
	clust(class configFile* conf, class imf_class* gimf, class iso_class_main* iso);
	~clust();
	//methods
	std::vector<float> generate_stars(double MAX_CLUSTER_MASS, double START_RND, double RND_RATIO, gsl_rng* random_number_gsl, std::vector<float> &BINARIES);
	void write_output(mat CLUST, double CURRENT_CLUSTER_MASS, int cur_ssp_index);
	double set_lower_mass_limit(class configFile* conf, int cur_ssp_index);
	void set_stars_gen_limits(double LOWER_MASS, double& START_RND, double& RND_RATIO, double& fCUR_CLUST_MASS);
	fcolvec optimal_generation_of_stars(float MAX_CLUSTER_MASS, double& CURRENT_CLUSTER_MASS, int cur_ssp_index);
	double get_m_next(double m,IMF *imf);

	//variables
	//random number generator
	gsl_rng*	random_number_gsl;
	gsl_rng*	gaussian_number_gsl;
	const gsl_rng_type * gsl_type;
	//pflamm-alteburg imf struct
	IMF imf;
	fcolvec SSP_Mmax;
	float SSP_Mmax_applied; // !!! dangerous

	//clases
	class configFile* conf;
	class imf_class* gimf;
	class iso_class_main* iso;
};
