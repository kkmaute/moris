#include "st_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace mtk
    {

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Pair::Mesh_Pair(
                Interpolation_Mesh* aInterpolationMesh,
                Integration_Mesh*   aIntegrationMesh,
                bool                aIsOwned)
                : mInterpolationMesh(aInterpolationMesh)
                , mIntegrationMesh(aIntegrationMesh)
                , mIsOwned(aIsOwned)
        {
            MORIS_ERROR(aInterpolationMesh, "Interpolation mesh must exist in a mesh pair.");
            /* FIXME should also check the integration mesh, but currently this breaks something in the mapper when it
            is forced to create a paired integration mesh for a given HMR interpolation mesh */

            // MORIS_ERROR(aIntegrationMesh, "Integration mesh must exist in a mesh pair.");
        }

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Pair::Mesh_Pair(const Mesh_Pair& aMeshPair)
                : mInterpolationMesh(aMeshPair.mInterpolationMesh)
                , mIntegrationMesh(aMeshPair.mIntegrationMesh)
                , mIsOwned(false)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Mesh_Pair::~Mesh_Pair()
        {
            if (mIsOwned)
            {
                delete mInterpolationMesh;
                delete mIntegrationMesh;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Interpolation_Mesh* Mesh_Pair::get_interpolation_mesh() const
        {
            return mInterpolationMesh;
        }

        //--------------------------------------------------------------------------------------------------------------

        Integration_Mesh* Mesh_Pair::get_integration_mesh() const
        {
            MORIS_ASSERT(mIntegrationMesh, "Integration mesh does not exist.");
            return mIntegrationMesh;
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}