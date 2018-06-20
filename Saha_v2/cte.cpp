#include "cte.hpp"

#include <string>

/****   PDG 2014  ****/
const double long hbar=6.58211928e-16; // eV.s
const double long k_B=8.6173324e-5; // Boltzmann Cte eV/K
const double long m_e=5.10998928e3; // mass of an electron eV/c²
const double long v_light=299792458; // speed of light m/s
const double long Ry=13.6; // Rydberg Cte eV
const double long PI=3.14159265359; //


/****  NIST  ****/
// Fe
const double long E_Fe_FeI = 7.9024681 ; // eV
const double long E_FeI_FeII = 16.19920 ;  // eV
const double long E_FeII_FeIII = 30.651 ; // eV

//Bi
const double long E_Bi_BiI = 7.9024681 ; // eV
const double long E_BiI_BiII = 16.19920 ;  // eV
const double long E_BiII_BiIII = 30.651 ; // eV


/*** Z(X) Fits - Gray p.514 ***/
/*** the first number is the order of the polyfit ***/
// Fe
const double long Z_o5_FeI[] = {
    5,
    -2.73036859, 17.35740822, -42.37110286,
    49.57143794, -28.01764044, 7.678
};

const double long Z_o5_FeII[]= {
    5,
    -0.56690705, 3.6394595, -9.00571824,
    10.75911276, -6.39959266, 3.22
};

const double long Z_o5_FeIII[]= {
    5,
    -2.73036859, 17.35740822, -42.37110286,
    49.57143794, -28.01764044, 7.678
};

extern const double long Z_o8_FeI[]= {
    8,
    2.97425285, -29.20706043, 122.60486877,
    -287.261923, 411.17795777, -369.25260338,
    204.72798848, -65.33775091, 11.019
};

extern const double long Z_o8_FeII[]= {
    8,
    0.58613126, -5.74927237, 24.12508298,
    -56.57303156, 81.22721035, -73.48678592,
    41.40827224, -13.78135319, 3.8815
};

extern const double long Z_o8_FeIII[]= {
    8,
    2.97425285, -29.20706043, 122.60486877,
    -287.261923, 411.17795777, -369.25260338,
    204.72798848,-65.33775091, 11.019
};


// Bi
const double long Z_o5_BiI[] = {
    5
};

const double long Z_o5_BiII[]= {
    5
};

const double long Z_o5_BiIII[]= {
    5
};

extern const double long Z_o8_BiI[]= {
    8
};

extern const double long Z_o8_BiII[]= {
    8
};

extern const double long Z_o8_BiIII[]= {
    8
};



int Z_number=0;
std::string Z_name;
bool verbose=false;
const char SEPARATOR=' ';
const std::string extname="dat"; // should be removed
unsigned int fit_order=0;
