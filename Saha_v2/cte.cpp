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
const double long E_Bi_BiI = 7.285516; // eV
const double long E_BiI_BiII = 16.703;  // eV 1986
const double long E_BiII_BiIII = 25.563; // eV 1958 and 1936


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
    5,  
    -0.76402244, 4.87474505, -11.98476836,
    14.2112398, -8.2170803, 2.5192
};

const double long Z_o5_BiII[]= {
    5, 
    -0.00761218, 0.17143794, -0.93468094, 
    2.21110868, -2.47385583, 1.08206667
};

const double long Z_o5_BiIII[]= {
    5, 
    -0.30729167, 1.94223485, -4.71983537,
    5.52338287, -3.15368228, 1.01586667
};

extern const double long Z_o8_BiI[]= {
    8, 
    0.8719308, -8.51186829, 35.5066636,
    -82.65078636, 117.55503217, -105.03542766, 
    58.1929232, -18.76221687, 3.4612
};

extern const double long Z_o8_BiII[]= {
    8, 
    0.23251488, -1.48718342, 2.94436785, 
    0.24466401, -8.70774612, 12.2333188,
    -6.15789227, -0.09628792, 0.8484
};

extern const double long Z_o8_BiIII[]= {
    8, 0.44080946, -4.21695371, 17.20904182,
    -39.12284263, 54.2380036, -47.09927963,
    25.24188494, -7.82445974, 1.4297
};



int Z_number=0;
std::string Z_name;
bool verbose=false;
const char SEPARATOR=' ';
const std::string extname="dat"; // should be removed
unsigned int fit_order=0;
