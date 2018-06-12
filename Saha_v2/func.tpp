#include "func.hpp"
#include "io.hpp"
#include "cte.hpp"

#include <iostream>
#include <cmath>
#include <vector>


/******* Electronic Pressure P_e ********/

template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::e_pressure(std::vector<_T> &N_e,std::vector<_T> &T) {
	double long dim=N_e.size();
	std::vector<_T> P_e;
	if (dim==T.size()) {
		for(unsigned int i=0;i<dim;i++) {
			P_e.push_back((N_e[i]*k_B*T[i]));
		}
	}
	else {
		std::cerr << "e_pressure(): vectors dim do not match !\n";
	}
	return P_e;
}


/******* Saha-Boltzmann Equation *******/

template<class _T, class _T1, class _T2>
std::vector<_T2> func<_T,_T1,_T2>::saha_Ne(
			std::vector<_T1> &Z1,
			std::vector<_T1> &Z0,
			std::vector<_T2> &Ne,
			std::vector<_T2> &T,
			int nb_e) {

	std::vector<_T2> pop_ratio;

	unsigned int dim=Ne.size();

	if (nb_e<1) {
		std::cerr << "saha(): invalid number of electron  !\n";
		nb_e=1;
	}
	if (dim!=Z0.size() || dim!=Z1.size() || dim!=T.size()) {
		std::cerr << "saha(): vector sizes do not match !\n";
	}
	else {
		//_T2 C1=pow(2.0*PI*pow(hbar,2.0)*pow(v_light,2.0)/(m_e*k_B),-1.5);  // raw formula
        _T2 C2= -1.5 * ( log10(2.0*PI) + 2.0*log10(hbar) + 2.0*log10(v_light) - log10(m_e)) ;

		for(unsigned int i=0;i<dim;i++) {
			_T2 z=-1.0; // to detect an error

			/* raw
            if (nb_e==24 && Ne[i]!=0) z = 2.0 * C1 * (Z1[i]/Z0[i]) * pow( T[i], 1.5)  * exp(-E_FeII_FeIII/(k_B * T[i])) / Ne[i];
            if (nb_e==25 && Ne[i]!=0) z = 2.0 * C1 * (Z1[i]/Z0[i]) * pow( T[i], 1.5)  * exp(-E_FeI_FeII/(k_B * T[i])) / Ne[i];
            */

            if (nb_e==25 && Ne[i]!=0) z = log10 (2.0) + C2 + Z1[i] - Z0[i] - log10(Ne[i]) + 1.5 * log10(k_B * T[i]) - (E_Fe_FeI+E_FeI_FeII)/(log(10.0) * k_B * T[i]);
            if (nb_e==24 && Ne[i]!=0) z = log10 (2.0) + C2 + Z1[i] - Z0[i] - log10(Ne[i]) + 1.5 * log10(k_B * T[i]) - (E_Fe_FeI+E_FeI_FeII+E_FeII_FeIII)/(log(10.0) * k_B * T[i]);

         pop_ratio.push_back(z);
		}
	}
    return pop_ratio;
}

template<class _T, class _T1, class _T2>
std::vector<_T2> func<_T,_T1,_T2>::saha_Pe(
			std::vector<_T1> &Z1,
			std::vector<_T1> &Z0,
			std::vector<_T2> &Pe,
			std::vector<_T2> &T,
			int nb_e) {

	std::vector<_T2> pop_ratio;

	unsigned int dim=Pe.size();

	if (nb_e<1) {
		std::cerr << "saha(): invalid number of electron  !\n";
		nb_e=1;
	}
	if (dim!=Z0.size() || dim!=Z1.size() || dim!=T.size()) {
		std::cerr << "saha(): vector sizes do not match !\n";
	}
	else {
        _T2 C2= -1.5 * ( log10(2.0*PI) + 2.0*log10(hbar) + 2.0*log10(v_light) - log10(m_e) ) ;

		for(unsigned int i=0;i<dim;i++) {
			_T2 z=-1.0; // to detect an error

            if (nb_e==25 && Pe[i]!=0) z = log10 (2.0) + C2 + Z1[i] - Z0[i] - log10(Pe[i]) + 2.5 * log10( k_B * T[i]) - (E_Fe_FeI+E_FeI_FeII)/(log(10.0) * k_B * T[i]);
            if (nb_e==24 && Pe[i]!=0) z = log10 (2.0) + C2 + Z1[i] - Z0[i] - log10(Pe[i]) + 2.5 * log10( k_B * T[i]) - (E_Fe_FeI+E_FeI_FeII+E_FeII_FeIII)/(log(10.0) * k_B * T[i]);

//             if (nb_e==25 && Pe[i]!=0) z = log10 (0.665) + Z1[i] - Z0[i] - log10(Pe[i]) + 2.5 * log10(  T[i]) - E_FeI_FeII/(log(10.0) * k_B * T[i]);
//             if (nb_e==24 && Pe[i]!=0) z = log10 (0.665) + Z1[i] - Z0[i] - log10(Pe[i]) + 2.5 * log10(  T[i]) - E_FeII_FeIII/(log(10.0) * k_B * T[i]);

         pop_ratio.push_back(z);
		}
	}
    return pop_ratio;
}


template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::Z_func(std::vector<_T> &T, int nb_e) {

    bool status=true;
    std::vector<_T1> Z;

   if (fit_order==8) {
        if (nb_e==26 && Z_o8_FeI[0]!=8) {std::cerr << "Z_Fe_8(): bad partition function for feI !\n"; status=false;}
        if (nb_e==25 && Z_o8_FeII[0]!=8) {std::cerr << "Z_Fe_8(): bad partition function for feII !\n"; status=false;}
        if (nb_e==24 && Z_o8_FeIII[0]!=8) {std::cerr << "Z_Fe_8(): bad partition function for feIII !\n"; status=false;}
    }
    if (fit_order==5) {
        if (nb_e==26 && Z_o5_FeI[0]!=5) {std::cerr << "Z_Fe_5(): bad partition function for feI !\n"; status=false;}
        if (nb_e==25 && Z_o5_FeII[0]!=5) {std::cerr << "Z_Fe_5(): bad partition function for feII !\n"; status=false;}
        if (nb_e==24 && Z_o5_FeIII[0]!=5) {std::cerr << "Z_Fe_5(): bad partition function for feIII !\n"; status=false;}
    }

    if (status) {

        if (verbose) std::cout << "Z_func(): computing partition function " << fit_order << ".\n";

        _T1 C=1.0/(log(10.0)*k_B); // 5040

        for(unsigned int i=0;i<T.size();i++) {
            _T1 zi=0.0;
            for(int p=fit_order;p!=-1;p--) { // higher powers come first
                if (fit_order==8) {
                    if (nb_e==26) zi+=Z_o8_FeI[fit_order+1-p]*pow(C/T[i],p);
                    if (nb_e==25) zi+=Z_o8_FeII[fit_order+1-p]*pow(C/T[i],p);
                    if (nb_e==24) zi+=Z_o8_FeIII[fit_order+1-p]*pow(C/T[i],p);
                }
                if (fit_order==5) {
                    if (nb_e==26) zi+=Z_o5_FeI[fit_order+1-p]*pow(C/T[i],p);
                    if (nb_e==25) zi+=Z_o5_FeII[fit_order+1-p]*pow(C/T[i],p);
                    if (nb_e==24) zi+=Z_o5_FeIII[fit_order+1-p]*pow(C/T[i],p);
                }
            }
            Z.push_back(zi);
        }
    }
    return Z;
}


template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::ratio_2elem(std::vector<_T> &pop_ratio) {
    std::vector<_T> res;
    for(unsigned int i=0;i<pop_ratio.size();i++) {
        _T x = pow(10.0,pop_ratio[i]);
         res.push_back(x/(1.0+x));
    }
    return res;
}


template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::ratio_3elem ( std::vector<_T> &pop_ratio1, std::vector<_T> &pop_ratio2 ) {
    std::vector<_T> res;
    unsigned int dim=pop_ratio1.size();
    if (dim!=pop_ratio2.size()) {
        std::cerr << "saha(): vector sizes do not match !\n";
    }
    else {
        for(unsigned int i=0;i<dim;i++) {
            _T x1 = pow(10.0,pop_ratio1[i]);
            _T x2 = pow(10.0,pop_ratio2[i]);
            res.push_back(x1/(1.0+x1+x2));
        }
    }
    return res;
}

template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::ratio_4elem ( std::vector<_T> &pop_ratio1, std::vector<_T> &pop_ratio2, std::vector<_T> &pop_ratio3 ) {
    std::vector<_T> res;
    unsigned int dim=pop_ratio1.size();
    if (dim!=pop_ratio2.size() || dim!=pop_ratio3.size()) {
        std::cerr << "saha(): vector sizes do not match !\n";
    }
    else {
        for(unsigned int i=0;i<dim;i++) {
            _T x1 = pow(10.0,pop_ratio1[i]);
            _T x2 = pow(10.0,pop_ratio2[i]);
			_T x3 = pow(10.0,pop_ratio3[i]);
            res.push_back(x1/(1.0+x1+x2+x3));
        }
    }
    return res;
}

template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::log10_vec(const std::vector<_T> &array) {
    std::vector<_T> res;
    for(unsigned int i=0;i<array.size();i++)  {
        _T x=array[i];
        if (x!=0) res.push_back(log10(x));
        else {
            res.push_back(0.0);
        }
    }
    return res;
}

template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::pow10_vec( const std::vector<_T> &array ) {
	std::vector<_T> res;
    for(unsigned int i=0;i<array.size();i++)  {
        _T x=array[i];
        res.push_back(pow(10.0,x));
    }
    return res;
}

template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::add_vector(const std::vector<_T> &X, const std::vector<_T> &Y) {
    std::vector<_T> res;
    if (X.size()>Y.size()) {
        for(unsigned int i=0;i<X.size();i++)
            res.push_back(static_cast<_T>(X[i])+static_cast<_T>(Y[i]));
    }
    else {
        for(unsigned int i=0;i<Y.size();i++)
            res.push_back(static_cast<_T>(X[i])+static_cast<_T>(Y[i]));
    }
    return res;
}

template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::lin_transform_vector ( const std::vector<_T> &X, _T1 a, _T2 b) {
	std::vector<_T> res;
        for(unsigned int i=0;i<X.size();i++)
            res.push_back(static_cast<_T>(a*X[i]+b));
    return res;
}


template<class _T, class _T1, class _T2>
std::vector<_T> func<_T,_T1,_T2>::inv_vector ( const std::vector<_T> &X) {
	std::vector<_T> res;
	for(unsigned int i=0;i<X.size();i++) {
		_T x=X[i];
		if (x!=0) res.push_back(1.0/x);
		else { res.push_back(0.0); }
	}
	return res;
}
