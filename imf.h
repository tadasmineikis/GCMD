#ifndef FIRST_TIME_IMF_H
#define FIRST_TIME_IMF_H

#include <string>
using namespace std;


class imf_class
{
	public:
		imf_class(int rows, int cols, float** array);
		~imf_class();

		double GenStar(double RandNum);
		double GenStar_MASS(double star_mass);
		double Cumulative_MASS(double star_mass);
		double Cumulative_NUM(double star_mass);
		
		int num_of_mass_intervals;
		float **MassSlope;

	private:
		//Methods
		//----------------------------
		bool Allocate();
		void Normalize();

		//Variables
		//----------------------------
		//imf describing parameters
		// starting mass point and slope
		
		double *join_constants_N, *integrated_imf_N, *join_constants_M, *integrated_imf_M;
		double *imf_gen_star_routine_N, *imf_gen_star_routine_2_N, *imf_gen_star_routine_3_N;
		double *imf_gen_star_routine_M, *imf_gen_star_routine_2_M, *imf_gen_star_routine_3_M;
		// mass ranges
		double Sum_N, Sum_M;
};

#endif
