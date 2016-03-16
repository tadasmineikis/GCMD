#include <iostream>
#include <cstdlib>
#include <fstream>
#include <float.h>
#include <vector>
#include <armadillo>
#define FIRST_TIME_FILTERS_N_H

#include "iso.h"
#include "lib.h"
#include "fileio.h"
#include "filters.h"

#define YY_col 17
#define PAD_col 21
#define GEN_col 21
#define PD_col

#define LINE_LENGTH 256

#define VERBOSE 0

using namespace std;
using namespace arma;

//----------------------------------------
// Main class-container
iso_class_main::iso_class_main(int num_of_files, string* file_names)
{
    // read number of input files
	i_iso = num_of_files;
	// create auxilary class'es and feed with file names
	iso = new iso_class*[i_iso];
	for(int i=0; i<i_iso; i++)
	{
		iso[i] = new iso_class(file_names[i]);
	}
}

iso_class_main::~iso_class_main()
{
    for(int i=0; i<i_iso; i++)
        delete iso[i];
    delete[] iso;
}
int iso_class_main::find_nearest_iso(float z)
{
    float min_d = 1;                //'distance' in z-untis, initialize to *large* value
    int min_iso_index = -1;         // index of 'right' isochrone
    for (int i=0; i<i_iso; i++)
    {
        if(abs(iso[i]->Z - z) < min_d)  // take modulus of z 'distance'
        {
            min_d = abs(iso[i]->Z - z);
            min_iso_index = i;
        }
    }
    if (min_iso_index == -1)
    {
        cerr << "Isochrone bank error? z=" << z << endl;
        for (int i=0; i<i_iso; i++)
            cerr <<"Z=" << iso[i]->Z << endl;
        cerr << "aborting..." << endl;
        exit(EXIT_FAILURE);
    }
    return min_iso_index;
}

float iso_class_main::find_limiting_lower_mass(float ssp_age, float ssp_z, float lim_value, int lim_filter)
{
    int nearest_z = find_nearest_iso(ssp_z); //index of nearest isochrone
    return iso[nearest_z]->find_lower_mass(lim_value, lim_filter, ssp_age);
}

float iso_class_main::find_limiting_upper_mass(float ssp_age, float ssp_z)
{
    int nearest_z = find_nearest_iso(ssp_z); //index of nearest isochrone
    return iso[nearest_z]->find_upper_mass(ssp_age);
}

//---------------------------------------
// auxilary isochrone class
iso_class::iso_class(string fname)
{
	//open isochrone file, check if ok
	//-----------------------------
	fstream in;
	in.open(fname.c_str(), ios::in);
	if(!in.is_open())
	{
	    cerr<< "File not found. [iso_class::iso_class(string fname" + fname + ")]" << endl;
	    exit(EXIT_FAILURE);
	}

	read_generic_isochrone(in);
	print();
	in.close();
}

iso_class::~iso_class()
{
    free(mass_range),free(pop_age);
}

void iso_class::print()
{
#if VERBOSE
	cout <<"ages:" << ages <<endl;
	cout <<"Z:" << Z <<endl;
	cout <<"Y:" << Y <<endl;
	cout <<"Feh:" << FeH <<endl;
	cout <<"aFe:" << aFe <<endl;
	cout << pop_data[ ages-1 ][ mass_range[ages-1] - 1] [ 0 ] << endl;
	cout << endl;
#endif 
}

double iso_class::find_lower_mass(float lim_value, int lim_filter, float log_ssp_age)
{
    int nearest_age = find_nearest_age(log_ssp_age);
	if(pop_data[nearest_age](0,lim_filter) < lim_value)
	{
		cerr << "Lower mass cut off lower than lowest available!" <<endl;
		cerr << "Using lowest mass available in isochrones instead!" << endl;
		cerr << "Old value: " << pow(10, lim_value) << "\t" << lim_value << "\t";
		lim_value = pop_data[nearest_age](0,enum_mass_i);
		cerr << "new value: " << pow(10, lim_value) << "\t" << lim_value << endl;
		return lim_value;
	}
	else
	{
    	for(int i=1; i<mass_range[nearest_age]; i++)
    	{
        	if(lim_value >= pop_data[nearest_age](i,lim_filter))
        	{
            	return LinInterpP(
                    pop_data[nearest_age](i - 1,0), pop_data[nearest_age](i,0),
                    pop_data[nearest_age](i - 1,lim_filter),	pop_data[nearest_age](i,lim_filter),
					lim_value);
        	}
    	}
	}
    return 200;
}

float iso_class::find_upper_mass(float log_ssp_age)
{
    int nearest_age = find_nearest_age(log_ssp_age);
    return pop_data[nearest_age](mass_range[nearest_age] -1,0);
}

int iso_class::find_nearest_age(float log_c_age)
{
    float min_d = 100;      //'distance' in age-untis, initialize to *large* value
    int min_iso_index = -1;         // index of 'right' isochrone
    for (int i=0; i<ages; i++)
    {
        // for skipping unessesary iterations
        if(abs(pop_age[i] - log_c_age) < min_d)  // take modulus of age 'distance'
        {
            min_d = abs(pop_age[i] - log_c_age);
            min_iso_index = i;
        }
    }
    if (min_iso_index == -1)
    {
        cerr << "Isochrone error? age=" << log_c_age << endl;
        for (int i=0; i<ages; i++)
            cerr <<"age=" << pop_age[i] << endl;
        cerr << "aborting..." << endl;
        exit(EXIT_FAILURE);
    }
    return min_iso_index;
}

rowvec iso_class::only_mass_interpolate( float log_mass_s, int age_index, int& last_mas_index)
{
	bool ERROR_FLAG = true;

    //interpolating mass
    rowvec result(num_of_columns + 2);
    if(log_mass_s < pop_data[age_index](0,0))
	{
	    //result = (float*) calloc(num_of_columns + 2, sizeof(float));
	    //stars outside lower boundary
	    //get values of last mass grid point
	    result.cols(0,num_of_columns-1) = pop_data[age_index].row(0);
	    result(0) = log_mass_s;

        result(num_of_columns) = Z;
		result(num_of_columns+1) = pop_age[age_index];
		ERROR_FLAG = false;
		return result;
    }
    if(log_mass_s > pop_data[age_index](mass_range[age_index] - 1,0))
	{
	    //result = (float*) calloc(num_of_columns + 2, sizeof(float));
	    //stars outside lower boundary
	    //get values of last mass grid point
	    //set all elements to zeros
	    result.cols(0,num_of_columns-1).fill(-512);
	    result(0) = log_mass_s;

        result(num_of_columns) = Z;
		result(num_of_columns+1) = pop_age[age_index];
		ERROR_FLAG = false;
		return result;
    }
    int i;
    for(i = last_mas_index; i<mass_range[age_index]; i++)
	{
		if(log_mass_s < pop_data[age_index](i,0))
		{
			//result = (float*) calloc(num_of_columns + 2, sizeof(float));
			
			//(float y2, float y1, float x2, float x1, float x)
			//y1 + (y2 - y1)/(x2 - x1)*(x - x1);
			result.cols(0,num_of_columns-1) = (pop_data[age_index].row(i)- pop_data[age_index].row(i-1))/(pop_data[age_index](i,0) - pop_data[age_index](i-1,0))*(log_mass_s - pop_data[age_index](i-1,0)) + pop_data[age_index].row(i-1);
			result(0)=log_mass_s;

			// log10(pop_data[age_index][i][0]),	log10(pop_data[age_index][i -1][0]),
			result(num_of_columns) = Z;
			result(num_of_columns+1) = pop_age[age_index];
			ERROR_FLAG = false;
			return result;
		}
	}
	last_mas_index = i;
	if(ERROR_FLAG)
	{
        if(log_mass_s < pop_data[age_index](mass_range[age_index] - 1,0) + FLT_EPSILON)
		{
			//result = (float*) calloc(num_of_columns + 2, sizeof(float));
			result.cols(0,num_of_columns-1) = pop_data[age_index].row(mass_range[age_index] - 1);
			result(0) = log_mass_s;
			
			// log10(pop_data[age_index][i][0]),	log10(pop_data[age_index][i -1][0]),
			result(num_of_columns) = Z;
			result(num_of_columns+1) = pop_age[age_index];
			ERROR_FLAG = false;
			return result;
		}

	}

	if(ERROR_FLAG  && log_mass_s > pop_data[age_index](mass_range[age_index] - 1,0))
	{
	    cerr << "Error: stellar mass to high for such age\n";
	    cerr << log_mass_s << "\t" << pop_data[age_index](mass_range[age_index] - 1,0) << endl;
	    return NULL;
	}

    cerr << "Unexpexted behavior, aborting..." << endl;
    exit(EXIT_FAILURE);
}

// reading generic isochrone type
void iso_class::read_generic_isochrone(std::iostream& in)
{
    // read number of ages
    in >> ages;

    //allocation of data storage
    pop_age = new float[ages];
    mass_range = new int[ages];
    pop_data = new mat[ages];

    num_of_columns = 0;
    for(int i=0; !in.eof()&& i<ages; i++)
    {
        char trash;
        in >> trash;
        if(in.eof())
            break;

        in >> mass_range[i] >> Z >> pop_age[i] >> num_of_columns;

        //take log values
        pop_age[i] = log10(pop_age[i]);
        logZ = log10(Z);

        pop_data[i]=mat(mass_range[i],num_of_columns);

        for(int j=0; j<mass_range[i]; j++)
            for(int k=0; k<num_of_columns; k++)
            {
                in >> pop_data[i](j,k);
            }
    }
}


//======================================================================
//          ERROR class funtions
//
//======================================================================
//error::error(string* files, int num_of_files_)
error::error(float c1, float c2)
{
    //int cols=4, SkipRows=1;
    //num_of_files = num_of_files_;
    //data = new DataFile(num_of_files, files, cols, SkipRows);
    CONST_1 = c1;
    CONST_2 = c2;

}

int error::find_nearest_error(float mag)
{
    int i;
    for(i=1; i<num_of_files; i++)
    {
        if(mag >= error_values[i])
            break;
    }

    if(i == num_of_files)
        return i - 1;
    else
        return i;
}

float error::generate_error(float mag)
{
    //float sigma = mag;//generate_sigma(mag);
    float temp = 0;//box_muller(mag, sigma);
    //return box_muller(mag, sigma);
    return temp;
}

void error::init_error_values(float* error_values_)
{
    error_values = error_values_;
    error_values_ = NULL;
}

float error::generate_sigma(float mag)
{

    return pow(10, (CONST_1*mag + CONST_2));
}

/*float error::generate_sigma(float mag, float rand_num)
{
    int nearest_e = find_nearest_error(mag);

    for(int i=1; i<data->Rows[nearest_e]; i++)
    {
        if(rand_num <= data->data_storage[nearest_e][i][3] )
        {
            float temp = LinInterp(
                    data->data_storage[nearest_e][i][0], data->data_storage[nearest_e][i - 1][0],
                    data->data_storage[nearest_e][i][3],	data->data_storage[nearest_e][i - 1][3],
					rand_num);
            return temp;
        }
    }

    return -999;
}
*/
