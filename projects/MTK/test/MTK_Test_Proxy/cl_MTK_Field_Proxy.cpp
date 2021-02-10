
#include "../MTK_Test_Proxy/cl_MTK_Field_Proxy.hpp"

#include "fn_dot.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {
        Field_Proxy::Field_Proxy(
                std::shared_ptr<mtk::Mesh_Manager>   aMeshManager,
                uint const                         & aMeshIndex,
                uint const                         & aDiscretizationMeshIndex )
        : Field( aMeshManager,
                aMeshIndex,
                aDiscretizationMeshIndex )
        {
        }

        // ----------------------------------------------------------------------------------------------

        Field_Proxy::~Field_Proxy()
        {

        }

        // ----------------------------------------------------------------------------------------------

    }
}

