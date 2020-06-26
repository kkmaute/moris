#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        Mesh_Manager::Mesh_Manager()
        {
        }

        //-------------------------------------------------------------------------

        Mesh_Manager::~Mesh_Manager()
        {
            uint tCounter = 0;

            for( bool tIsOwned : mIsOwned )
            {
                if( tIsOwned )
                {
                    MORIS_ERROR( mInterpolationMesh(tCounter) != nullptr, "Interpolation mesh already deleted." );

                    delete mInterpolationMesh(tCounter);

                    MORIS_ERROR( mIntegrationMesh(tCounter) != nullptr, "Integration mesh already deleted." );

                    delete mIntegrationMesh(tCounter);
                }
                tCounter++;
            }

            mInterpolationMesh.clear();
            mIntegrationMesh.clear();
        }

        uint
        Mesh_Manager::register_mesh_pair(
                Interpolation_Mesh* aInterpMesh,
                Integration_Mesh*   aIntegrationMesh,
                bool                aIsOwned)
        {
            mInterpolationMesh.push_back(aInterpMesh);
            mIntegrationMesh.push_back(aIntegrationMesh);

            // set flag which tells if mesh is owned by manager. If it is owned it ha to be deleted
            mIsOwned.push_back(aIsOwned);

            return mInterpolationMesh.size()-1;
        }

        void
        Mesh_Manager::get_mesh_pair(
                moris::moris_index     aPairIndex,
                Interpolation_Mesh * & aInterpMesh,
                Integration_Mesh   * & aIntegrationMesh)
        {
            MORIS_ASSERT(
                    aPairIndex < (moris::moris_index) mInterpolationMesh.size() &&
                    aPairIndex < (moris::moris_index) mIntegrationMesh.size(),
                    "Mesh_Manager::get_integration_mesh: requested mesh does not exist.");

            aInterpMesh      = mInterpolationMesh(aPairIndex);
            aIntegrationMesh = mIntegrationMesh(aPairIndex);
        }

        Interpolation_Mesh*
        Mesh_Manager::get_interpolation_mesh(moris::moris_index aMeshIndex)
        {
            MORIS_ASSERT( aMeshIndex < (moris::moris_index) mInterpolationMesh.size(),
                    "Mesh_Manager::get_integration_mesh: requested mesh does not exist.");

            return mInterpolationMesh(aMeshIndex);
        }

        Integration_Mesh*
        Mesh_Manager::get_integration_mesh(moris::moris_index aMeshIndex)
        {
            MORIS_ASSERT( aMeshIndex < (moris::moris_index) mIntegrationMesh.size(),
                    "Mesh_Manager::get_integration_mesh: requested mesh does not exist.");

            return mIntegrationMesh(aMeshIndex);
        }
    }
}

