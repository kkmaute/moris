#ifndef MORIS_ST_GEN_PROPERTY_PARAMETERS_HPP
#define MORIS_ST_GEN_PROPERTY_PARAMETERS_HPP

#include "st_GEN_Field_Parameters.hpp"
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * This struct contains additional parameters that are used by properties.
         */
        struct Property_Parameters
        {
            PDV_Type mPDVType = PDV_Type::UNDEFINED;  //! The type of PDV that this property will be assigned to
            bool mInterpolationPDV = true;            //! If the PDV is defined on the interpolation mesh (always true for now)
            Cell<std::string> mPDVMeshSetNames = {};  //! Mesh set names for assigning PDVs
            Matrix<DDUMat> mPDVMeshSetIndices = {{}}; //! Mesh set indices for assigning PDVs
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
