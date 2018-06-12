#ifndef CTE_HPP
#define CTE_HPP

#include <string>

/*** Constantes ***/

extern const double long hbar;
extern const double long k_B;
extern const double long m_e;
extern const double long v_light;
extern const double long Ry;
extern const double long PI;


/*** E(X) ***/

extern const double long E_Fe_FeI;
extern const double long E_FeI_FeII;
extern const double long E_FeII_FeIII;


/*** Z(X) Fitted ***/
extern const double long Z_o5_FeI[];
extern const double long Z_o5_FeII[];
extern const double long Z_o5_FeIII[];

extern const double long Z_o8_FeI[];
extern const double long Z_o8_FeII[];
extern const double long Z_o8_FeIII[];


/*** g(X) ***/

extern const double long g_Fe;
extern const double long g_FeI;
extern const double long g_FeII;
extern const double long g_FeIII;


extern bool verbose;
extern const char SEPARATOR;
extern const std::string extname;
extern unsigned int fit_order;

#endif
