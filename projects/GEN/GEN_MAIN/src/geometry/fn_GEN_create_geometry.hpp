//
// Created by christopherson on 4/13/20.
//

#ifndef MORIS_FN_CREATE_GEOMETRY_HPP
#define MORIS_FN_CREATE_GEOMETRY_HPP

#endif //MORIS_FN_CREATE_GEOMETRY_HPP

#include "cl_GEN_Geometry_Analytic.hpp"
#include "cl_Param_List.hpp"
#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Creates an instance of the specified Geometry class and returns a shared pointer to it
         *
         * @param aGeometryParameterList Parameter list for creating a geometry class
         * @param aADVs Reference to the initial adv vector
         * @param aLibrary pointer to library for loading user-defined functions
         * @return Pointer to specific Geometry class
         */
        std::shared_ptr<Geometry_Analytic> create_geometry(ParameterList aGeometryParameterList, Matrix<DDRMat>& aADVs, std::shared_ptr<moris::Library_IO> aLibrary = nullptr);

    }
}
