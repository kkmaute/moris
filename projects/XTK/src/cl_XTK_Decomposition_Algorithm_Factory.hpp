/*
 * cl_XTK_Decomposition_Algorithm_Factory.hpp  
 * 
 *  Created on: Oct  01, 2021 
 *      Author: Keenan Doble
 */
#ifndef SRC_cl_XTK_Decomposition_Algorithm_Factory
#define SRC_cl_XTK_Decomposition_Algorithm_Factory

#include "cl_XTK_Enums.hpp"
#include <memory>

namespace xtk
{
    class Decomposition_Algorithm;

    // free function
    std::shared_ptr<Decomposition_Algorithm>
    create_decomposition_algorithm(enum Subdivision_Method aSubdivisionMethod);

}



#endif /* cl_XTK_Decomposition_Algorithm_Factory.hpp */