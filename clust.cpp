#include "armadillo"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <float.h>
#include <iostream>
#include <iomanip>

#include "clust.h"
#include "imf.h"
#include "configFile.h"
#include "iso.h"
#define FIRST_TIME_FILTERS_N_H
#include "filters.h"


//pflamm-alteburg dependencies
#include "imf_pflamm.h"

using namespace arma;

clust::clust(configFile* conf_t, imf_class* gimf_t, iso_class_main* iso_t)
{
	//copy classes
	conf=conf_t;
	gimf=gimf_t;
	iso=iso_t;

	//if seed was not set, set it from clock
	if ( !conf->RANDOM_GEN_SEED_SET )
	{
		conf->RANDOM_GEN_SEED = time(NULL);
	}

    //gsl_rng_env_setup();
	if ( conf->sampling_method.at("stochastic") )
	{
		cerr<<"Sampling method: stochastic"<<endl;
	}
	else if ( conf->sampling_method.at("optimal") )
	{
		gsl_type =gsl_rng_mt19937;
		cerr<<"Sampling method: optimal"<<endl;
		gaussian_number_gsl = gsl_rng_alloc(gsl_type);
		gsl_rng_set(gaussian_number_gsl, conf->RANDOM_GEN_SEED);

		char *imf_string;
		imf_string = (char *) malloc(1001*sizeof(char));
		stringstream imf_params;
		//write imf defining parameters into string
		for (int i=0; i<gimf->num_of_mass_intervals;i++)
		{
			imf_params<<gimf->MassSlope[i][0]<<"(pow:"<<gimf->MassSlope[i][2]<<")";
		}
		imf_params<<gimf->MassSlope[gimf->num_of_mass_intervals-1][1];
		strcpy(imf_string,imf_params.str().c_str());
		//sent to pfalmm-alteburg routine of imf initialisation
		imf_init_user(&imf,imf_string);

		double mmaxphys = gimf->MassSlope[gimf->num_of_mass_intervals-1][1];
		SSP_Mmax=zeros<fcolvec>(conf->ssp_num_of_ssp);
		for(int ii=0; ii<conf->ssp_num_of_ssp; ii++ )
		{
			SSP_Mmax(ii)=imf_norm_cl_pawk11_Mmax(&imf,pow(10,conf->ssp_log_mass[ii]),mmaxphys);
		}
		free(imf_string);
	}
	cerr << "SEED:" << conf->RANDOM_GEN_SEED << endl;
}

clust::~clust()
{
	//release random number generator
	gsl_rng_free (random_number_gsl);
}

double clust::set_lower_mass_limit(configFile* conf, int cur_ssp_index)
{
	double lower_mass_cut = 1e9;
	if(conf->filter_cut != 0)
	{
	    lower_mass_cut = iso->find_limiting_lower_mass(conf->ssp_log_age[cur_ssp_index], conf->ssp_z[cur_ssp_index], conf->value_cut, conf->filter_cut) + FLT_EPSILON;
		lower_mass_cut = pow(10, lower_mass_cut);
		//cerr<<conf->lower_mass_cut << "\t" << conf->upper_mass_cut << endl;
	}
	//for cuts in mass
	else
	{
		lower_mass_cut = conf->value_cut;
	}
	return lower_mass_cut;
}

void clust::set_stars_gen_limits(double LOWER_MASS, double& START_RND, double& RND_RATIO, double& fCUR_CLUST_MASS)
{
	double RND_MAX=1;
	if(RND_MAX -  conf->ssp_integrate_until > 0 && conf->start_integration_flag)
    {
		LOWER_MASS = gimf->GenStar( RND_MAX -  conf->ssp_integrate_until );
		START_RND = RND_MAX -  conf->ssp_integrate_until;
		RND_RATIO = conf->ssp_integrate_until;
	}
	else if (conf->start_integration_flag )
	{
		LOWER_MASS = gimf->GenStar( 0 );
		START_RND = 0;
		RND_RATIO = RND_MAX;
	}
	else if ( !conf->start_integration_flag )
	{
		START_RND = gimf->Cumulative_NUM(LOWER_MASS);
		RND_RATIO = RND_MAX - START_RND;
	}
    fCUR_CLUST_MASS =  gimf->Cumulative_MASS( LOWER_MASS ) ;
}

std::vector<float> clust::generate_stars(double MAX_CLUSTER_MASS, double START_RND, double RND_RATIO, gsl_rng* rnd_gsl, std::vector<float> &vBINARIES)
{
	const int RESERVE_NUM_OF_STARS = 1024*64;
	//fcolvec STARS(RESERVE_NUM_OF_STARS);

	std::vector<float> vSTARS;
	//std::vector<float> vBINARIES;
	vSTARS.reserve(RESERVE_NUM_OF_STARS);
	vBINARIES.reserve(RESERVE_NUM_OF_STARS);

    double CURRENT_CLUSTER_MASS=0;
    for(int i=0; CURRENT_CLUSTER_MASS < MAX_CLUSTER_MASS; ++i)
    {
		vSTARS.push_back( gimf->GenStar( START_RND + RND_RATIO*gsl_rng_uniform (rnd_gsl) ) );
		CURRENT_CLUSTER_MASS += vSTARS[i];
		if(gsl_rng_uniform (rnd_gsl)< conf->BINARY_FRACTION)
		{
			vBINARIES.push_back( vSTARS[i]*gsl_rng_uniform (rnd_gsl) );
			CURRENT_CLUSTER_MASS += vBINARIES[i];
		}
		else
		{
			vBINARIES.push_back(1e10);
		}
    }
    return vSTARS;
}

fcolvec clust::optimal_generation_of_stars(float MAX_CLUSTER_MASS, double& CURRENT_CLUSTER_MASS, int cur_ssp_index)
{
	fcolvec STARS(10000000);
	double Mecl,mmaxphys,m,M_diced;

	int i;

	Mecl = MAX_CLUSTER_MASS;

	float gauss_rnd_num=gsl_ran_gaussian(gaussian_number_gsl,1);
	mmaxphys = gimf->MassSlope[gimf->num_of_mass_intervals-1][1];
	float tmp_Mmax = pow(10,log10(SSP_Mmax(cur_ssp_index))+conf->optimal_sampling_sigma*gauss_rnd_num);
	SSP_Mmax_applied= tmp_Mmax >= mmaxphys?
			mmaxphys : tmp_Mmax;

	imf_norm_cl_pawk11_fast(&imf,Mecl,mmaxphys,SSP_Mmax_applied);

	m = imf.m_max;
	M_diced = 0;
	i = 1;
	int STARS_INDEX=0;
	while(M_diced <= Mecl){
		M_diced += m;
		STARS(STARS_INDEX)=m;
		STARS_INDEX++;
		//fprintf(stdout,"%d %g %g\n",i,m,M_diced);
		if(imf_int_mxi(&imf,imf.m[0],m) < imf.m[0])
			break;
		m = get_m_next(m,&imf);
		i++;
	}
	CURRENT_CLUSTER_MASS=M_diced;
	STARS.resize(STARS_INDEX);
	return  sort(log10(STARS));
}

double clust::get_m_next(double m,IMF *imf){
  double a,b,c,mb;
  const double err = 1.e-8;

  a = imf->m[0];
  c = m;
  b = (c+a)/2;
  mb = imf_int_mxi(imf,b,m);
  while(fabs((mb-b)/mb) > err)
    {
      mb = imf_int_mxi(imf,b,m);
      if(mb < b)
	c = b;
      else
	a = b;
      b = (c+a)/2;
    }
  return b;
}

void clust::write_output(mat CLUST, double CURRENT_CLUSTER_MASS, int cur_ssp_index)
{
	rowvec FILTERS(enum_LAST);

			uvec good_index=arma::find(CLUST.cols(enum_mass_c, enum_mass_c) < -500,1,"first");
			//cerr<<"end"<<endl;
			mat CLUST_LUM;
			if(good_index.n_rows > 1)
			{
				//cerr<<"start0"<<endl;
				//cerr<<good_index(0) <<" "<<  CLUST(0,enum_mass_c) << " " <<  CLUST(0,enum_mass_i) << endl;
				CLUST_LUM = CLUST.submat(0,enum_mass_c+1,good_index(0)-1,enum_LAST-1);
				//cerr<<"end0"<<endl;
			}
			else
			{
				//cerr<<"start1"<<endl;
				CLUST_LUM = CLUST.submat(0,enum_mass_c+1,CLUST.n_rows-1,enum_LAST-1);
				//cerr<<"end1"<<endl;
			}
			FILTERS.cols(enum_mass_c+1, enum_LAST-1) = sum(exp10(-0.4*CLUST_LUM),0);
			CLUST_LUM.clear();
			FILTERS.cols(enum_mass_i,enum_mass_i)=sum(exp10(CLUST.cols(enum_mass_i, enum_mass_i)),0);
			FILTERS.cols(enum_mass_c,enum_mass_c)=sum(exp10(CLUST.cols(enum_mass_c, enum_mass_c)),0);

			cout << setprecision(4); 	//decimal digits
			cout << log10(CURRENT_CLUSTER_MASS) << " " << FILTERS(enum_mass_i) << " ";
			if(1. - conf->ssp_integrate_until > 0)
				cout << gimf->GenStar( 1.- conf->ssp_integrate_until ) << " " ;
			else
				cout << gimf->GenStar( 0 ) << " " ;
			if(conf->sampling_method.at("optimal"))
			{
				cout << SSP_Mmax(cur_ssp_index) << " "<< SSP_Mmax_applied << " ";
			}
			else if(conf->sampling_method.at("stochastic"))
			{
				cout << SSP_Mmax_applied << " ";
			}
			else
			{
				cerr <<"Error!"<< endl;
				exit(EXIT_FAILURE);
			}
			cout << conf->ssp_integrate_until <<" ";
			cout << 1. - gimf->Cumulative_MASS(gimf->GenStar(gimf->Cumulative_NUM(conf->upper_mass_cut) - conf->ssp_integrate_until )) << " ";
			cout << CLUST(0,iso->iso[0]->num_of_columns + 1) << " " << CLUST(0,iso->iso[0]->num_of_columns) << " ";

			for(int enum_iterator=enum_mass_c+1; enum_iterator<enum_LAST; ++enum_iterator)
				cout << -2.5*log10(FILTERS(enum_iterator)/CURRENT_CLUSTER_MASS) <<" ";
			cout <<endl;

	CLUST.clear();
	FILTERS.clear();
}
