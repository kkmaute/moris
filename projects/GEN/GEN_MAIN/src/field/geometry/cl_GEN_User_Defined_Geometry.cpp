#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        User_Defined_Geometry::User_Defined_Geometry(
                Matrix<DDRMat>            aConstants,
                Field_Function            aFieldFunction,
                Geometry_Field_Parameters aParameters)
                : Field(aConstants, aParameters)
                , User_Defined_Field(aConstants, aFieldFunction, aParameters)
                , Geometry(aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
