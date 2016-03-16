#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <float.h>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <omp.h>

#include "configFile.h"
#include "lib.h"
#include "filters.h"

configFile::configFile(char *f_name)
{
	set_init_values();
	//opening file procedures
	fstream input_file;
	input_file.open(f_name, ios::in);

	if(!input_file.is_open())
	{
		cerr<<"Open operation for file: " << f_name <<" generated error. Aborting..." << endl;
		exit(EXIT_FAILURE);
	}

	// Config file rules:
	// 1. Lines starting with # are ignored

	string input_str;

	for(;;)
	{
		input_file >> input_str;

		if(input_str.size() == 0 || input_file.eof())
		{
			break;
		}
		// check if line contents are comments, otherwise interpret
		else if(input_str.find_first_of("#") < input_str.npos )
            getline(input_file, input_str);
        else if(input_str.find_first_of("#") == input_str.npos )
		{
			if(input_str.find("iso", 0) != input_str.npos)
			{
			   	input_file >> iso_num_of_files;
			   	iso_files = new string[iso_num_of_files];
			   	for (int i=0; i<iso_num_of_files; i++)
			   	{
			   	    input_file >> iso_files[i];
			   	}
            }
			else if(input_str.find("out_fits", 0) != input_str.npos)
			{
			   	input_file >> FITS_SIZE_X >> FITS_SIZE_Y >> FITS_FILE_NAME;
			   	OUT_FITS=true;
			   	init_out_fits();
            }
			else if(input_str.find("imf", 0) != input_str.npos)
			{
				string define_sampling_method;
				input_file >> define_sampling_method >> optimal_sampling_sigma;
			    input_file >> imf_num_of_intervals >> imf_num_of_cols;
			    imf_def_array = (float**)matrix(imf_num_of_intervals, imf_num_of_cols, sizeof(float));

			    for(int i=0; i<imf_num_of_intervals; i++)
                    for(int j=0; j<imf_num_of_cols; j++)
                        input_file >> imf_def_array[i][j];

                //initialize dictionary of method names
                int list_index=2;
                string list_sampling_methods[]={
					"stochastic",
					"optimal"
				};
				for(int i=0; i<list_index; i++)
					sampling_method[list_sampling_methods[i]];
                sampling_method.at(define_sampling_method) = true;
            }
            else if(input_str.find("limit", 0) != input_str.npos)
			{
			    input_file >> filter_cut >> value_cut;
            }
            else if(input_str.find("ssp_integration", 0) != input_str.npos)
			{
				start_integration_flag = true;
			    input_file >> ssp_integrate_until;
            }
            else if(input_str.find("threads", 0) != input_str.npos)
			{
				int threads;
				input_file >> threads ;
			    omp_set_num_threads(threads);
            }
            else if(input_str.find("ssp", 0) != input_str.npos)
			{
			    input_file >> ssp_num_of_ssp;

			    ssp_log_age = new float[ssp_num_of_ssp];
			    ssp_log_mass = new float[ssp_num_of_ssp];
			    ssp_z = new float[ssp_num_of_ssp];
                ssp_copies = new int[ssp_num_of_ssp];

			    for(int i=0; i<ssp_num_of_ssp; i++)
			    {
			        input_file >> ssp_log_age[i] >>  ssp_z[i] >> ssp_log_mass[i] >> ssp_dist_mod >> ssp_copies[i];
			    }
			    //output header into std
			    if (sampling_method["optimal"])
			    	    cout <<"logMass_i_Total Mass_i MinStarMass OptimalMaxStarMass MaxStarMass nFrac mFrac logAge Z";
				else if(sampling_method["stochastic"])
						cout <<"logMass_i_Total Mass_i MinStarMass MaxStarMass nFrac mFrac logAge Z";
				for(int i=enum_mass_c+1; i<enum_LAST; i++)
					cout <<" "<<enum_filter_names[i];
				cout << endl;

            }
            else if(input_str.find("galemo_results", 0) != input_str.npos)
			{
				galemo_ssp_flag = true;
				string gin_file_name;
			    input_file >> ssp_num_of_ssp>>gin_file_name;
			    ssp_dist_mod = 0;

			    ssp_log_age = new float[ssp_num_of_ssp];
			    ssp_age = new int[ssp_num_of_ssp];
			    ssp_log_mass = new float[ssp_num_of_ssp];
			    ssp_z = new float[ssp_num_of_ssp];
                ssp_radius = new float[ssp_num_of_ssp];
			    ssp_azimut = new float[ssp_num_of_ssp];
			    ssp_x = new float[ssp_num_of_ssp];
			    ssp_y = new float[ssp_num_of_ssp];
                ssp_copies = new int[ssp_num_of_ssp];


                gsl_rng*	random_number_gsl;
				const gsl_rng_type * gsl_type = gsl_rng_mt19937;

                random_number_gsl = gsl_rng_alloc (gsl_type);
				gsl_rng_set(random_number_gsl, time(NULL));

				fstream fgin;
				fgin.open(gin_file_name.c_str(),ios::in);
				if (fgin.is_open() == false)
				{
					cerr<<"Failure to open file: "<< gin_file_name.c_str() << endl;
					exit(EXIT_FAILURE);
				}
				//read header
				getline(fgin, input_str);
				cerr<<"READING input ssp file..." <<endl;
			    for(int i=0; i<ssp_num_of_ssp; ++i)
			    {
			        fgin >> ssp_radius[i] >> ssp_azimut[i] >> ssp_log_mass[i] >> ssp_age[i] >> ssp_z[i];
			        ssp_copies[i] = 1;
			        if (ssp_age[i] == 0)
			        {
						ssp_log_age[i] = 6.6+gsl_rng_uniform(random_number_gsl)*0.4;
					}
					else
					{
						ssp_log_age[i] = log10(ssp_age[i]+gsl_rng_uniform(random_number_gsl)-0.5)+6;
					}
			    }
			   /* for(int i=0; i<ssp_num_of_ssp; ++i)
			    {
					cerr << ssp_radius[i] <<" " << ssp_azimut[i]<<" " << ssp_log_mass[i]<<" "<< ssp_age[i]<<"/ " << ssp_log_age[i] << " " << ssp_z[i] << endl;
					cerr<< ssp_num_of_ssp << endl;
				}*/
				cerr<<"FINISHED READING input ssp file..." <<endl;
			    fgin.close();

            }
            else if(input_str.find("out", 0) != input_str.npos)
			{
			    input_file >> out_file_name;
				out_cmd = true;
			    out_file.open(out_file_name.c_str(), ios::out);

			    //write header for cmd file
				if(galemo_ssp_flag)
				{
					out_file <<"#r a x y age a_age z a_z mass_i_p mass_i_s";
					for(int i=enum_mass_c+1; i<enum_LAST; i++)
						out_file <<" p_"<<enum_filter_names[i];
					for(int i=enum_mass_c+1; i<enum_LAST; i++)
						out_file <<" s_"<<enum_filter_names[i];
					for(int i=enum_mass_c+1; i<enum_LAST; i++)
						out_file <<" o_"<<enum_filter_names[i];
					out_file <<" ssp_index"<<endl;
				}
				else
				{
					out_file <<"#age z mass_i mass_c";
					for(int i=enum_mass_c+1; i<enum_LAST; i++)
						out_file <<" "<<enum_filter_names[i];
					out_file <<" ssp_index"<<endl;
				}
            }
            else if(input_str.find("error", 0) != input_str.npos)
			{
                input_file >> error_c1 >> error_c2;
            }
            else if(input_str.find("BINARY_FRACTION", 0) != input_str.npos)
			{
                input_file >> BINARY_FRACTION;
                cerr << BINARY_FRACTION << endl;
            }
            else if(input_str.find("seed", 0) != input_str.npos)
			{
                input_file >> RANDOM_GEN_SEED;
                RANDOM_GEN_SEED_SET = true;
            }
			else if(input_str.find("predefined", 0) != input_str.npos)
			{
				string predefined_file_name;
                input_file >> predefined_file_name >> PREDEFINED_NUM_OF_LINES >> PREDEFINED_Z >> PREDEFINED_LOG_AGE;
                fstream predefined_file;
                predefined_file.open(predefined_file_name.c_str(),ios::in);
                string predefined_line;
                char comment;
                comment=predefined_file.peek();
                while (comment == '#')
                {
					predefined_file.ignore( 1024, '\n');
					comment=predefined_file.peek();
				}
				PREDEFINED_LIST = (float*) calloc(PREDEFINED_NUM_OF_LINES, sizeof(float));
				for (int i=0; i < PREDEFINED_NUM_OF_LINES; i++)
						predefined_file >> PREDEFINED_LIST[i];

				predefined_list_flag = true;

            }
            else if(input_str.find("perlin-dust", 0) != input_str.npos)
			{
			    input_file >> ssp_num_of_ssp;

			    ssp_log_age = new float[ssp_num_of_ssp];
			    ssp_log_mass = new float[ssp_num_of_ssp];
			    ssp_z = new float[ssp_num_of_ssp];
                ssp_radius = new float[ssp_num_of_ssp];
			    ssp_azimut = new float[ssp_num_of_ssp];
                ssp_gas_mass = new float[ssp_num_of_ssp];
			    //ssp_dist_mod = new float[ssp_num_of_ssp];

			    for(int i=0; i<ssp_num_of_ssp; i++)
			    {
			        input_file >> ssp_log_age[i] >>  ssp_z[i] >> ssp_log_mass[i] >> ssp_gas_mass[i] >> ssp_dist_mod;
			        ssp_radius[i] = 0, ssp_azimut[i] = 0;
			    }

            }
		}
	}
	input_file.close();
}

configFile::~configFile()
{
	if(OUT_FITS)
		write_out_fits();
    delete[] iso_files;
    //free_matrix((void**)imf_def_array);

    //ssp
    delete[] ssp_log_age;
    delete[] ssp_log_mass;
    delete[] ssp_z;
    delete[] ssp_radius;
    delete[] ssp_azimut;
    //delete[] ssp_dist_mod;

    //out
    out_file.close();

}

void configFile::set_init_values()
{
	galemo_ssp_flag = false;
	start_integration_flag = false;
}

void configFile::output_cmd(arma::mat CLUST, arma::mat BIN_CLUST, int ssp_index, gsl_rng* rnd_gsl, int binaries)
{
	//cerr<<CLUST.n_cols<< " "<< enum_LAST+1 <<" "<< index <<endl;
	//cerr<<CLUST(0,enum_LAST+1)<<endl;
	//CLUST.save("dump.dat",arma::raw_ascii);
	static int cmd_id = 0;

	stringstream cmd_str;

	float cell_w, x, y, rad_noise, azi_noise;
	float total_mag;

	if(galemo_ssp_flag)
	{
		//cerr<<"start"<<endl;
		//cerr<< cmd_id << endl;
		for(unsigned int i=0; i<CLUST.n_rows; ++i)
		{
			total_mag=-999;
			if(CLUST(i,enum_mass_c) != -512 && BIN_CLUST(i,enum_mass_c)== -512)
			{
				total_mag = CLUST(i,filter_cut);
			}
			else if(CLUST(i,enum_mass_c) != -512 && BIN_CLUST(i,enum_mass_c) != -512 )
				total_mag = -2.5*log10(pow(10,-0.4*CLUST(i,filter_cut))+pow(10,-0.4*BIN_CLUST(i,filter_cut)));
			else if(CLUST(i,enum_mass_c) == -512 && BIN_CLUST(i,enum_mass_c) != -512 )
				total_mag = BIN_CLUST(i,filter_cut);
			//if(CLUST(i,enum_mass_c) != -512 && CLUST(i,filter_cut) < value_cut)
			if( (CLUST(i,enum_mass_c) != -512 || BIN_CLUST(i,enum_mass_c) != -512) && (total_mag < value_cut) )
			{
				if(ssp_radius[ssp_index]==0)
					cell_w=2*arma::datum::pi;
				else
					cell_w=2*arma::datum::pi/(ssp_radius[ssp_index]/0.05*6);
				rad_noise = 0.05*gsl_rng_uniform(rnd_gsl)-0.025;
				azi_noise = cell_w*(gsl_rng_uniform(rnd_gsl)-0.5);

				x=(ssp_radius[ssp_index]+rad_noise)*cos(ssp_azimut[ssp_index]+azi_noise);
				y=(ssp_radius[ssp_index]+rad_noise)*sin(ssp_azimut[ssp_index]+azi_noise);

				cmd_str << setprecision(3);
				cmd_str << ssp_radius[ssp_index];        //RADIUS
				cmd_str << " " << ssp_azimut[ssp_index]; // AZIMUTH
				cmd_str << " " << x << " " << y;         // X, Y
				cmd_str << setprecision(4);
				cmd_str << " " << ssp_log_age[ssp_index] << " " << CLUST(i,enum_LAST+1); // AGE, ACTUAL ISOCHRONE AGE
				cmd_str << " " << ssp_z[ssp_index] << " " << CLUST(i,enum_LAST);         // METALLICITY, ACTUAL ISOCHRONE METALLICITY
				cmd_str << setprecision(4);
				cmd_str << " " << CLUST(i,enum_mass_i) <<" "<< BIN_CLUST(i,enum_mass_i); // PRIMARY MASS, SECONDARY
				if(CLUST(i,enum_mass_c) != -512){
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter)
					{
						cmd_str  <<  " " <<CLUST(i,i_filter);                                //FILTERS OF PRIMARY STAR
					}
				}
				else{
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter)
					{
						cmd_str  <<   " \"\" ";                                  //FILTERS
					}
				}
				if(BIN_CLUST(i,enum_mass_c) != -512){
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter) {
						cmd_str  <<  " " <<BIN_CLUST(i,i_filter);                                //FILTERS OF SECONDARY STAR
					}
				}
				else{
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter) {
						cmd_str  <<  " \"\" ";                                //FILTERS
					}
				}
				if(BIN_CLUST(i,enum_mass_c) != -512 && CLUST(i,enum_mass_c) != -512){
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter) {
						cmd_str  <<  " " <<-2.5*log10(pow(10,-0.4*CLUST(i,i_filter))+pow(10,-0.4*BIN_CLUST(i,i_filter))); //COMBINED FILTERS
					}
				}
				else if ( (CLUST(i,enum_mass_c) == -512) && (BIN_CLUST(i,enum_mass_c) != -512) ){
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter) {
						cmd_str  <<  " " <<BIN_CLUST(i,i_filter);
					}
				}
				else if ( (CLUST(i,enum_mass_c) != -512) && (BIN_CLUST(i,enum_mass_c) == -512) ){
					for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; ++i_filter) {
						cmd_str  <<  " " <<CLUST(i,i_filter);
					}
				}
				cmd_str <<" "<< cmd_id <<"\n";                                           // SSP ID
			}
		}
		//cerr<<"finish"<<endl;
	}
	else
	{
		cmd_str << setprecision(5);
		for(unsigned int i=0; i<CLUST.n_rows; i++)
		{
			cmd_str << ssp_log_age[ssp_index];
			cmd_str << " " << ssp_z[ssp_index];
			cmd_str << " " << CLUST(i,enum_mass_i) << " " << CLUST(i,enum_mass_c);
			for (int i_filter=enum_mass_c+1; i_filter<enum_LAST; i_filter+=1)
			{
				cmd_str  <<  " " <<CLUST(i,i_filter);
			}
			cmd_str <<" "<< cmd_id <<"\n";
		}
	}
	out_file << cmd_str.str();
	cmd_id+=1;
}


void configFile::print_input(char* fn)
{
}

void configFile::init_out_fits()
{
	FITS_ARRAY = (float**) matrix(FITS_SIZE_X,FITS_SIZE_Y,sizeof(float));
	//memset(FITS_ARRAY[0], 0,sizeof(float)*FITS_SIZE_X*FITS_SIZE_Y);
	MISSING_FLUX=0;
}

void configFile::out_fits(float** cluster_stars_properties, int filter_index, int index)
{
	int X, Y;
	for (int i=0; i<index; i++)
	{
		X=int((FITS_SIZE_X)*((float)rand()/RAND_MAX-DBL_EPSILON));
		Y=int((FITS_SIZE_Y)*((float)rand()/RAND_MAX-DBL_EPSILON));
		FITS_ARRAY[X][Y]+= pow(10,-0.4*cluster_stars_properties[i][filter_index]);
	}
	for (int i=0; i<FITS_SIZE_X; i++)
		for(int j=0; j<FITS_SIZE_Y; j++)
			FITS_ARRAY[i][j]+=1;
}

void configFile::write_out_fits()
{
	/*long naxes[2] = { FITS_SIZE_X, FITS_SIZE_Y };
	long nelements = naxes[0]*naxes[1];
	int status = 0;
	fitsfile* FITS_FILE;
	//remove previous file
	remove(FITS_FILE_NAME.c_str());
	//create new file
	string header("Z_FLUX");
//	char MISSING_FLUX_HEADER[]={"missing flux"};
	if(!fits_create_file(&FITS_FILE, FITS_FILE_NAME.c_str(), &status))
	{
		if(!fits_create_img(FITS_FILE, FLOAT_IMG, 2, naxes, &status))
		{
			if(!fits_write_img(FITS_FILE, TFLOAT, 1, nelements, FITS_ARRAY[0], &status))
			{
//				if (! fits_update_key(FITS_FILE, TFLOAT, header.c_str(), &MISSING_FLUX, &MISSING_FLUX_HEADER, &status) )
//				{
					if (! fits_close_file(FITS_FILE, &status) )
						status =0;
//				}
			}
		}
	}
	else
		cerr << "Error producing fits file." <<endl;*/


}
