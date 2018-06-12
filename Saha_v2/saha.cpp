#include <iostream>
#include <vector>
#include <string>
#include <ctime> //timer
#include <algorithm> //min max

#include "io.hpp"
#include "func.hpp"
#include "cte.hpp"
#include "py_wrapper.hpp"
//#include "plot.hpp" // unable to link SDL2 & SDL2_ttf


int main ( int argc, char** argv ) {
    clock_t start=clock();

// init class
    func<> _func; // assume double long
    io<> _io;
    py_wrapper _wrapper;

// parse cmd line and get params
    _io.parse_cmd ( argc, argv );

// init arrays
    std::vector<std::vector<double long> > data;

    std::vector<double long> Z_FeI, Z_FeII, Z_FeIII;
    std::vector<double long> res_FeI, res_FeII, res_FeIII, res_FeI_tmp0, res_FeI_tmp1; // mid computations

    std::vector<double long> ratio_FeII_I_II; // ne(FeII) / ( ne(FeI) + ne(FeII) )
    std::vector<double long> ratio_FeI_I_II; // ne(FeI) / ( ne(FeI) + ne(FeII) )
    std::vector<double long> ratio_FeI_I_II_III; // ne(FeI) / ( ne(FeI) + ne(FeII) + ne(Fe III) )
    std::vector<double long> ratio_FeII_I_II_III; // ne(FeII) / ( ne(FeI) + ne(FeII) + ne(Fe III) )
    std::vector<double long> ratio_FeIII_I_II_III; // ne(FeIII) / ( ne(FeI) + ne(FeII) + ne(Fe III) )

    std::vector<double long> opt_depth;


// on lit les données du modele fort.8
    _io.read_atmos ( data );

// on calcul les polynomes des fonctions de partition avec T
    Z_FeI=_func.Z_func ( data[1], 26 ); // 26 e-
    Z_FeII=_func.Z_func ( data[1], 25 ); // 25 e- e.g. Fe II
    Z_FeIII=_func.Z_func ( data[1], 24 ); // 24 e- e.g. Fe III...

// Compute pop ratios
    res_FeII=_func.saha_Ne ( Z_FeII, Z_FeI, data[3], data[1], 25 );
    res_FeIII=_func.saha_Ne ( Z_FeIII, Z_FeI, data[3], data[1], 24 );

    ratio_FeI_I_II=_func.inv_vector(_func.lin_transform_vector(_func.pow10_vec(res_FeII),1.0,1.0));

    ratio_FeII_I_II=_func.ratio_2elem ( res_FeII );
    ratio_FeII_I_II_III=_func.ratio_3elem ( res_FeII, res_FeIII );
    ratio_FeIII_I_II_III=_func.ratio_3elem ( res_FeIII, res_FeII );

    ratio_FeI_I_II_III=_func.add_vector(
                        _func.lin_transform_vector( ratio_FeII_I_II_III , -1.0, 1.0),
                        _func.lin_transform_vector( ratio_FeIII_I_II_III, -1.0, 0.0) );

    opt_depth=_func.log10_vec ( data[0] );

// write data
    io<> _io_1, _io_2, _io_3, _io_4,  _io_5;

    _io_1.output=_io.output;
    _io_2.output=_io.output;
    _io_3.output=_io.output;
    _io_4.output=_io.output;
    _io_5.output=_io.output;

    _io_1.write_csv ( opt_depth, ratio_FeI_I_II );
    _io_2.write_csv ( opt_depth, ratio_FeII_I_II );
    _io_3.write_csv ( opt_depth, ratio_FeI_I_II_III );
    _io_4.write_csv ( opt_depth, ratio_FeII_I_II_III );
    _io_5.write_csv ( opt_depth, ratio_FeIII_I_II_III );


// plot with py.matplotlib wrapper
    _wrapper.plot ( _io_1.output_stamp, "Population ratio as function of Optical depth", "$log(\\\\tau_{5000}$)","$\\\\frac{n_{FeI}}{n_{FeI}+n_{FeII}}$" );
    _wrapper.plot ( _io_2.output_stamp, "Population ratio as function of Optical depth", "$log(\\\\tau_{5000})$","$\\\\frac{n_{FeII}}{n_{FeI}+n_{FeII}}$" );
    _wrapper.plot ( _io_3.output_stamp, "Population ratio as function of Optical depth", "$log(\\\\tau_{5000}$)","$\\\\frac{n_{FeI}}{n_{FeI}+n_{FeII}+n_{FeIII}}$" );
    _wrapper.plot ( _io_4.output_stamp, "Population ratio as function of Optical depth", "$log(\\\\tau_{5000}$)","$\\\\frac{n_{FeII}}{n_{FeI}+n_{FeII}+n_{FeIII}}$" );
    _wrapper.plot ( _io_5.output_stamp, "Population ratio as function of Optical depth", "$log(\\\\tau_{5000}$)","$\\\\frac{n_{FeIII}}{n_{FeI}+n_{FeII}+n_{FeIII}}$" );

    std::cout << "elapsed time: " << ( clock()-start ) / ( ( float ) CLOCKS_PER_SEC ) << " s\n";
    return 0;
}
