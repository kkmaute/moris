#include "cl_MTK_Mesh_Manager.hpp"

namespace moris
{
    namespace mtk
    {
        //-------------------------------------------------------------------------
        
        Mesh_Manager::Mesh_Manager()
        {
        }

        //-------------------------------------------------------------------------

        Mesh_Manager::~Mesh_Manager()
        {
        }

        //-------------------------------------------------------------------------

        uint
        Mesh_Manager::register_mesh_pair(
                Interpolation_Mesh* aInterpolationMesh,
                Integration_Mesh*   aIntegrationMesh,
                bool                aIsOwned)
        {
            // Create new mesh pair
            Mesh_Pair tMeshPair;
            tMeshPair.mInterpolationMesh = aInterpolationMesh;
            tMeshPair.mIntegrationMesh = aIntegrationMesh;
            tMeshPair.mIsOwned = aIsOwned;
            
            // Push back new mesh
            mMeshPairs.push_back(tMeshPair);

            return mMeshPairs.size() - 1;
        }

        //-------------------------------------------------------------------------

        void
        Mesh_Manager::get_mesh_pair(
                moris_index          aPairIndex,
                Interpolation_Mesh * & aInterpolationMesh,
                Integration_Mesh   * & aIntegrationMesh)
        {
            MORIS_ASSERT(
                    aPairIndex < (moris_index) mMeshPairs.size(),
                    "Mesh_Manager::get_mesh_pair: requested mesh pair does not exist.");

            aInterpolationMesh = mMeshPairs(aPairIndex).mInterpolationMesh;
            aIntegrationMesh = mMeshPairs(aPairIndex).mIntegrationMesh;
        }

        //-------------------------------------------------------------------------

        Interpolation_Mesh*
        Mesh_Manager::get_interpolation_mesh(moris_index aMeshIndex)
        {
            MORIS_ASSERT( aMeshIndex < (moris_index) mMeshPairs.size(),
                    "Mesh_Manager::get_interpolation_mesh: requested mesh does not exist.");

            return mMeshPairs(aMeshIndex).mInterpolationMesh;
        }

        //-------------------------------------------------------------------------

        Integration_Mesh*
        Mesh_Manager::get_integration_mesh(moris_index aMeshIndex)
        {
            MORIS_ASSERT( aMeshIndex < (moris_index) mMeshPairs.size(),
                    "Mesh_Manager::get_integration_mesh: requested mesh does not exist.");

            return mMeshPairs(aMeshIndex).mIntegrationMesh;
        }

        //-------------------------------------------------------------------------
    }
}

