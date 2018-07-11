/*
 * cl_Function.cpp
 *
 *  Created on: Oct 28, 2016
 *      Author: doble
 */

#include "cl_Function.hpp"


moris::tools::Function::Function()
{

}


moris::real
moris::tools::Function::evaluate_function(moris::Mat<moris::real>  &aCoords)
{
    return this->compute_func_val(aCoords);
}

moris::Mat<moris::real>
moris::tools::Function::evaluate_dxdp(moris::Mat<moris::real>  &aCoords)
{
    return this->compute_func_dxdp(aCoords);
}

