#include <string>
#include <armadillo>
#include <map>
#include <gsl/gsl_rng.h>

using namespace std;
class configFile
{
public:

	//constructor
	configFile(char *f_name);
	~configFile();

	/*functions*/
	void print_input(char* fn);
	void set_init_values();
	void output_cmd(arma::mat, arma::mat, int ssp_index, gsl_rng* random_number_gsl, int binaries);
	//fits handling
	void init_out_fits();
	void out_fits(float** cluster_stars_properties, int filter_index, int index);
	void write_out_fits();


	//variables
	int imf_num_of_intervals;        //number defining of mass intervals
	int imf_num_of_cols;        //number of cols in read-line
	float**imf_def_array;       //imf parametru masyvas

    int iso_num_of_files;       //number of isochrones
    string* iso_files;          //isochrones file names container

    float value_cut;
    int filter_cut;
    double lower_mass_cut, upper_mass_cut;
    
    //binaries
    float BINARY_FRACTION;

    int ssp_num_of_ssp;
    float *ssp_log_age;
    int *ssp_age;
    float *ssp_log_mass;
    float *ssp_z;
    float ssp_dist_mod;
    float* ssp_gas_mass;
    float *ssp_radius;
    float *ssp_azimut;
    float *ssp_x, *ssp_y;
	int *ssp_copies;
	bool galemo_ssp_flag;

    string out_file_name;
    fstream out_file;
	bool out_cmd;
    char **out_buffer, *out_ptr;
    int out_buffer_size;

    //error tables
    int num_of_error_tables;
    string* error_files;
    float *error_values;
    float error_c1, error_c2;

    //perlin-dust parameters
    float Zeff; // dulkiu skales aukstis
    float deMagnitude; // amplitudes gesimo faktorius, t.y. sekanti iteracija turi 1/deMag amplitudes
    int nStep; //itearciju kiekis

	//ssp integration
	bool start_integration_flag;
    float ssp_integrate_until;

    //fits output
    int FITS_SIZE_X, FITS_SIZE_Y, CELL_WIDTH;
    string FITS_FILE_NAME;
	bool OUT_FITS;
	float** FITS_ARRAY;
	float MISSING_FLUX;

	//predefined
	bool predefined_list_flag;
	float* PREDEFINED_LIST;
	int PREDEFINED_NUM_OF_LINES;
	float PREDEFINED_Z, PREDEFINED_LOG_AGE;

	//random generator
	unsigned long int RANDOM_GEN_SEED;
	bool RANDOM_GEN_SEED_SET;

	//cluster stars sampling
	map <string,bool> sampling_method;
	float optimal_sampling_sigma;

};
