//
// Created by christopherson on 4/13/20.
//

#ifndef MORIS_FN_CREATE_GEOMETRY_HPP
#define MORIS_FN_CREATE_GEOMETRY_HPP

#endif //MORIS_FN_CREATE_GEOMETRY_HPP

#include "cl_GEN_Geometry.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Creates an instance of the specified Geometry class and returns a shared pointer to it
         *
         * @param aGeometryParameterList Parameter list for creating a geometry class
         * @param aADVs Reference to the initial adv vector
         * @return Pointer to specific Geometry class
         */
        std::shared_ptr<Geometry> create_geometry(ParameterList aGeometryParameterList, Matrix<DDRMat>& aADVs);
    }
}
