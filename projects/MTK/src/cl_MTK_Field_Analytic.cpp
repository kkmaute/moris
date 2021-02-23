
#include "cl_MTK_Field_Analytic.hpp"

#include "fn_dot.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {
        Field_Analytic::Field_Analytic(
                mtk::Mesh_Pair * aMeshPairs,
                uint const     & aDiscretizationMeshIndex )
        : Field( aMeshPairs,
                 aDiscretizationMeshIndex )
        {
        }

        // ----------------------------------------------------------------------------------------------

        Field_Analytic::~Field_Analytic()
        {

        }

        // ----------------------------------------------------------------------------------------------

    }
}

