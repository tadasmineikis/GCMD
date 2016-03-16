#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "lib.h"
#include "imf.h"

using namespace std;

#define IMF_COL 5
//---------------------------------------
// IMF class constructor
imf_class::imf_class(int rows, int cols, float** array)
{
	// take array adress
    MassSlope = array;
    num_of_mass_intervals = rows;
	//Number
    integrated_imf_N = new double[rows];
    join_constants_N = new double[rows];
    imf_gen_star_routine_N = new double[rows];
    imf_gen_star_routine_2_N = new double[rows];
    imf_gen_star_routine_3_N = new double[rows];
    //Mass
    integrated_imf_M = new double[rows];
    join_constants_M = new double[rows];
    imf_gen_star_routine_M = new double[rows];
    imf_gen_star_routine_2_M = new double[rows];
    imf_gen_star_routine_3_M = new double[rows];
    array = NULL;


    for(int i=0; i<num_of_mass_intervals; i++)
	{
		//Number
		join_constants_N[i] = pow(1./MassSlope[i][0], MassSlope[i][2]);
		for(int j= i-1; j>=0; j--)
		{
			join_constants_N[i] *= pow(MassSlope[j][1]/MassSlope[j][0], MassSlope[j][2]);
		}
		//Mass
		join_constants_M[i] = pow(1./MassSlope[i][0], MassSlope[i][2]);
		for(int j= i-1; j>=0; j--)
		{
			join_constants_M[i] *= pow(MassSlope[j][1]/MassSlope[j][0], MassSlope[j][2]);
		}
	}

	Normalize();
	//for(int i=0; i<num_of_mass_intervals; i++)
	//	cout <<MassSlope[i][0] << "\t"<<MassSlope[i][1]<<"\t" <<join_constants_N[i]/Sum_N <<endl;
	
	//cout << endl;
	
	//for(int i=0; i<num_of_mass_intervals; i++)
	//	cout <<MassSlope[i][0] << "\t"<<MassSlope[i][1]<<"\t" <<join_constants_M[i]/Sum_M <<endl;

}
//-------------------------------------------
// IMF class destructor
imf_class::~imf_class()
{
	free_matrix((void**)MassSlope);
	delete[] integrated_imf_N;
	delete[] join_constants_N;
	delete[] imf_gen_star_routine_N;
	delete[] imf_gen_star_routine_2_N;
	delete[] imf_gen_star_routine_3_N;
	
	delete[] integrated_imf_M;
	delete[] join_constants_M;
	delete[] imf_gen_star_routine_M;
	delete[] imf_gen_star_routine_2_M;
	delete[] imf_gen_star_routine_3_M;
}

void imf_class::Normalize()
{
	Sum_N = 0;
	for(int i=0; i< num_of_mass_intervals; i++)
	{
	    // integrating mass of segment
		integrated_imf_N[i] = (pow(MassSlope[i][1], ( MassSlope[i][2] + 1. )) -
				pow((MassSlope[i][0]), ( MassSlope[i][2] + 1.) ) )/ (MassSlope[i][2] + 1.) *
				join_constants_N[i];
		Sum_N += integrated_imf_N[i];
	}
	for(int i=0; i< num_of_mass_intervals; i++)
	{
		integrated_imf_N[i] /= Sum_N;
		imf_gen_star_routine_N[i] =  pow(MassSlope[i][0], (MassSlope[i][2] +1.));
		imf_gen_star_routine_2_N[i] = (MassSlope[i][2] +1.)*Sum_N/join_constants_N[i];
		imf_gen_star_routine_3_N[i] = 1./(MassSlope[i][2] +1.);

	}
	//Mass
	Sum_M = 0;
	for(int i=0; i< num_of_mass_intervals; i++)
	{
	    // integrating mass of segment
		integrated_imf_M[i] = (pow(MassSlope[i][1], ( MassSlope[i][2] + 2. )) -
				pow((MassSlope[i][0]), ( MassSlope[i][2] + 2.) ) )/ (MassSlope[i][2] + 2.) *
				join_constants_M[i];
		Sum_M += integrated_imf_M[i];
	}
	for(int i=0; i< num_of_mass_intervals; i++)
	{
		integrated_imf_M[i] /= Sum_M;
		imf_gen_star_routine_M[i] =  pow(MassSlope[i][0], (MassSlope[i][2] +2.));
	}
}
//-------------------------------------------
// Takes as input random number and generates
// star mass according to defined IMF
double imf_class::GenStar(double RandNum)
{
	double Mass = 0;
	double sub = 0;
	for(int i=0; i<num_of_mass_intervals; i++)
	{
		sub += integrated_imf_N[i];
		if( RandNum <= sub )
		{
			sub -= integrated_imf_N[i];
			Mass = imf_gen_star_routine_N[i] + (RandNum - sub)*imf_gen_star_routine_2_N[i];
			Mass = pow(Mass, imf_gen_star_routine_3_N[i]);
			break;
		}
	}
	if(RandNum >= 1)
	{
	    Mass = MassSlope[num_of_mass_intervals - 1][1];
	}
	if ( Mass == 0 )
	{
		cerr << "Error generating mass! [float imf_class::GenStar(float RandNum=" << RandNum << " " <<sub <<")]" <<endl;
		exit(EXIT_FAILURE);
	}
	return Mass;
}

double imf_class::Cumulative_MASS(double star_mass)
{
	double cumul_mass = 0;
	for (int i=0; i<num_of_mass_intervals; i++)
	{
		cumul_mass += integrated_imf_M[i];
		if(star_mass <= MassSlope[i][1])
		{
			cumul_mass -= integrated_imf_M[i];
			cumul_mass += (pow(star_mass, (MassSlope[i][2] +2.)) - imf_gen_star_routine_M[i])
			*join_constants_M[i]/Sum_M/(MassSlope[i][2] + 2.);
			break;
		}
	}
	return cumul_mass;
}

double imf_class::GenStar_MASS(double star_mass)
{
	for (int i=0; i<num_of_mass_intervals; i++)
	{
		if(star_mass <= MassSlope[i][1])
		{
			return pow(star_mass, (MassSlope[i][2] +1.))*join_constants_M[i]/Sum_M;
			break;
		}
	}
	return -1; //error
}

double imf_class::Cumulative_NUM(double star_mass)
{
	double cumul_NUM = 0;
	for (int i=0; i<num_of_mass_intervals; i++)
	{
		cumul_NUM += integrated_imf_N[i];
		if(star_mass <= MassSlope[i][1])
		{
			cumul_NUM -= integrated_imf_N[i];
			cumul_NUM += (pow(star_mass, (MassSlope[i][2] +1.)) - imf_gen_star_routine_N[i])
			*join_constants_N[i]/Sum_N/(MassSlope[i][2] + 1.);
			break;
		}
	}
	return cumul_NUM;
}
#undef IMF_COL
