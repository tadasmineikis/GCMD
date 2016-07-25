#define ACS 0
#define UBV 1

#ifndef FIRST_TIME_FILTERS_H
#define FIRST_TIME_FILTERS_H

#if ACS == 1
enum filters
{
	enum_mass_i,
	enum_mass_c,
	enum_F435W,
	enum_F555W,
	enum_F606W,
	enum_F625W,
	enum_F775W,
	enum_F814W,  	    	   	  
	enum_LAST
};
#elif UBV == 1
enum filters
{
	enum_mass_i,
	enum_mass_c,
	enum_U,
	enum_B,
	enum_V,
	enum_R,
	enum_I,
	enum_J,
	enum_H,
	enum_K,   	    	   	  
	enum_LAST
};
#endif
#endif

#ifndef FIRST_TIME_FILTERS_N_H
#define FIRST_TIME_FILTERS_N_H

#if ACS == 1
string enum_filter_names[] = 
{
	"mass_i",
	"mass",
	"F435W",
	"F555W",
	"F606W",
	"F625W",
	"F775W",
	"F814W"

};

#elif UBV == 1
string enum_filter_names[] = 
{
	"mass_i",
	"mass",
	"U",
	"B",
	"V",
	"R",
	"I",
	"J",
	"H",
	"K"

};
#endif
#endif
