/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * st_GEN_Property_Parameters.hpp
 *
 */

#ifndef MORIS_ST_GEN_PROPERTY_PARAMETERS_HPP
#define MORIS_ST_GEN_PROPERTY_PARAMETERS_HPP

#include "st_GEN_Field_Parameters.hpp"
#include "GEN_Data_Types.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * This struct contains additional parameters that are used by properties.
         */
        struct Property_Parameters
        {
            Cell<std::string> mDependencyNames = {};  //! Names of the dependencies of this property
            PDV_Type mPDVType = PDV_Type::UNDEFINED;  //! The type of PDV that this property will be assigned to
            bool mInterpolationPDV = true;            //! If the PDV is defined on the interpolation mesh (always true for now)
            Matrix<DDUMat> mPDVMeshSetIndices = {{}}; //! Mesh set indices for assigning PDVs
            Cell<std::string> mPDVMeshSetNames = {};  //! Mesh set names for assigning PDVs
        };

        /**
         * This is a struct used to simplify \ref moris::ge::Property constructors. It contains all property parameters.
         */
        struct Property_Field_Parameters : Field_Parameters, Property_Parameters
        {
        };
    }
}

#endif //MORIS_ST_GEN_PROPERTY_PARAMETERS_HPP

