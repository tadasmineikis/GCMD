#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <time.h>
#include <vector>
#include <string>
#include <float.h>
#include <armadillo>
#include <omp.h>

#include "iso.h"
#include "imf.h"
#include "configFile.h"
#include "clust.h"
#include "lib.h"

using namespace std;

template <typename T> vector<size_t> sort_indexes(const vector<T> &v);

int main(int argc, char *argv[])
{
	// initialise random number generator
	// check for config file
	if(argc < 2)
	{
	    cerr<<"Not all parameters suplied." <<endl;
	    exit(EXIT_FAILURE);
	}
	configFile* conf = new configFile(argv[1]);
	//input name for isochrone library
	iso_class_main* iso = new iso_class_main(conf->iso_num_of_files, conf->iso_files);
	//input name for imf definition file
	imf_class* m = new imf_class(conf->imf_num_of_intervals, conf->imf_num_of_cols, conf->imf_def_array);

	clust* cls = new clust(conf, m, iso);

	arma::wall_clock timer;
	timer.tic();


	const gsl_rng_type * gsl_type = gsl_rng_mt19937;
	gsl_rng* random_number_gsl = gsl_rng_alloc (gsl_type);
	gsl_rng_set(random_number_gsl, conf->RANDOM_GEN_SEED);

	cerr<<"ENETERING LOOP..."<<endl;
	for(int cur_ssp_index=0; cur_ssp_index<conf->ssp_num_of_ssp; ++cur_ssp_index)
	{
		for(int cur_ssp_copy=0; cur_ssp_copy<conf->ssp_copies[cur_ssp_index]; ++cur_ssp_copy)
		{
			cerr<<"CALLING lower mas cut..."<<endl;
			double lower_mass_cut = cls->set_lower_mass_limit(conf, cur_ssp_index);
			double fcur_cluster_mass, start_rnd, rnd_ratio;
			cls->set_stars_gen_limits(lower_mass_cut, start_rnd, rnd_ratio, fcur_cluster_mass);
			double mass_to_generate=pow(10,conf->ssp_log_mass[cur_ssp_index])*((double)1.-fcur_cluster_mass);
			std::vector<float> vSTARS;
			std::vector<float> vBINARIES;


			cerr<<"GENERATING stars..."<<endl;
			cerr<<"BINNARY MASS FRACTION: " << conf->BINARY_FRACTION <<endl;
			if(conf->sampling_method.at("stochastic"))
			{
				//cerr<<start_rnd << " " << rnd_ratio << endl;
				vSTARS = cls->generate_stars(mass_to_generate, start_rnd, rnd_ratio, random_number_gsl, vBINARIES);
				//SSP_Mmax_applied=STARS(STARS.n_rows-1);
			}
			else
			{
				cerr<<"Not implemented!"<<endl;
				exit(EXIT_FAILURE);
			}

			/*std::sort (vSTARS.begin(), vSTARS.end());
			if ( vSTARS[0] < lower_mass_cut*0.9 )
			{
				cerr<<"Something went wrong.. Writing file of stellar masses dump.dat" << endl;
				cerr<<lower_mass_cut << " " << vSTARS[0] << endl;
				fstream fdump;
				fdump.open("dump.dat",ios::out);
				for(unsigned int i=0; i<vSTARS.size(); ++i)
				{
					fdump<<vSTARS.at(i)<<endl;
				}
				fdump.close();
				//CLUST.save("dump.dat", raw_ascii);
				exit(EXIT_FAILURE);
			}*/
			int z_index = iso->find_nearest_iso(conf->ssp_z[cur_ssp_index]);
			int age_index = iso->iso[z_index]->find_nearest_age(conf->ssp_log_age[cur_ssp_index]);
			int last_mass_index = 0;
			mat CLUST = mat(vSTARS.size(), iso->iso[0]->num_of_columns+2);

			vector<size_t> idx_main= sort_indexes(vSTARS);
			for(unsigned int i=0; i<idx_main.size(); i++)
			{
				CLUST.row(idx_main[i]) = iso->iso[z_index]->only_mass_interpolate(log10(vSTARS[idx_main[i]]), age_index, last_mass_index);
			}
			idx_main.clear();
			
			/*cerr<<"binaries:"<< vBINARIES.size() <<endl;
			for(int i=0; i<vBINARIES.size(); i++){
				cerr<<vBINARIES[i]<< endl;
			}*/
			
			last_mass_index = 0;
			mat CLUST_BNR = mat(vBINARIES.size(), iso->iso[0]->num_of_columns+2);
			vector<size_t> idx_bnr= sort_indexes(vBINARIES);
			for(unsigned int i=0; i<idx_bnr.size(); i++)
			{
				CLUST_BNR.row(idx_bnr[i]) = iso->iso[z_index]->only_mass_interpolate(log10(vBINARIES[idx_bnr[i]]), age_index, last_mass_index);
			}
			idx_bnr.clear();
			//for(unsigned int i=0 ; i<vSTARS.size(); ++i)
			//{
			//	CLUST.row(i) = iso->iso[z_index]->only_mass_interpolate(log10(vSTARS[i]), age_index, last_mass_index);
			//}
			vSTARS.clear();
			vBINARIES.clear();
			conf->output_cmd(CLUST, CLUST_BNR, cur_ssp_index, random_number_gsl, 0);
			//conf->output_cmd( cur_ssp_index, random_number_gsl, 1);
			CLUST.clear();
			CLUST_BNR.clear();
		}
	}

	cerr << "took " << timer.toc() << " seconds" << endl;

    delete iso;
    delete m;
    delete conf;
    return EXIT_SUCCESS;
}

template <typename T> vector<size_t> sort_indexes(const vector<T> &v) 
{
  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
