#include <string>
#include <armadillo>
#include "imf.h"
using namespace std;
using namespace arma;

#define MAX_AGE 400

//=================================
//

//=================================
// Auxilary isochrone class
class iso_class
{
	public:
		iso_class(string fname);
		//iso_class(float c1, float c2);
		~iso_class();
	//methods
	//------------------
		void print();
		float* interpolate(double mass_s, float age_s);
		float* mean_age_interpolate( float log_mass_s, float age_s);
		rowvec only_mass_interpolate( float log_mass_s, int age_index, int& last_mass_index);
		int find_nearest_age(float c_age);
		double find_lower_mass(float lim_value, int lim_filter, float log_ssp_age);
		float find_upper_mass(float log_ssp_age);
		float integration_low_part(int age_index, float HIGH_MASS_CUT, int filter_index, class imf_class* m);

		float Z; 	//metal fraction by mass
		float logZ;
        int num_of_columns;

	private:
	//methods
	//------------------
        void read_generic_isochrone(std::iostream& in);
		void read_yy_par(std::iostream& in);
		void read_yy_data(std::iostream& in);
		void read_padova_par(std::iostream& in);
		void read_padova_data(std::iostream& in);

	//variables
	//------------------

		float Y;	//helium fraction by mass
		float FeH;	//[Fe/H]
		float aFe;	//[a/Fe]

		mat* pop_data; //isochrone
		float* pop_age;   //age of isochrone
		int* mass_range;
		int ages;
};

//=================================
// Main class-container
class iso_class_main
{
	public:
		iso_class_main(int num_of_files, string* file_names);
		~iso_class_main();

		iso_class** iso;
		float* interpolate(float mass_s, float age_s, float z_s);
		float* mean_z_isochrone(float mass_s, float age_s, float z_s);
		int find_nearest_iso(float z);
		float find_limiting_lower_mass(float ssp_age, float ssp_z, float lim_value, int lim_filter);
		float find_limiting_upper_mass(float ssp_age, float ssp_z);

		int i_iso;
};

class error
{
    public:
    //error(string* files, int num_of_files);
    error(float c1, float c2);
    ~error();

    float generate_error(float mag);
    float generate_sigma(float mag);
    int find_nearest_error(float mag);
    void init_error_values(float* error_values_);

    int num_of_files;
    float* error_values;
    class DataFile* data;
    float CONST_1, CONST_2;
};
