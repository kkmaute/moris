#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

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
            Mesh_Pair tMeshPair(aInterpolationMesh, aIntegrationMesh, false);
            
            // Push back new pair
            mMeshPairs.push_back(tMeshPair);

            // Ownership begins after adding
            mMeshPairs(mMeshPairs.size() - 1).mIsOwned = aIsOwned;

            return mMeshPairs.size() - 1;
        }

        //-------------------------------------------------------------------------

        uint Mesh_Manager::register_mesh_pair(Mesh_Pair& aMeshPair)
        {
            // Push back new pair
            mMeshPairs.push_back(aMeshPair);

            // Ensure correct ownership
            mMeshPairs(mMeshPairs.size() - 1).mIsOwned = aMeshPair.mIsOwned;
            aMeshPair.mIsOwned = false;

            return mMeshPairs.size() - 1;
        }

        //-------------------------------------------------------------------------

        const Mesh_Pair& Mesh_Manager::get_mesh_pair(moris_index aPairIndex)
        {
            MORIS_ASSERT( aPairIndex < (moris_index) mMeshPairs.size(),
                          "Mesh_Manager::get_mesh_pair: requested mesh pair does not exist.");
            return mMeshPairs(aPairIndex);
        }

        //-------------------------------------------------------------------------

        Mesh_Pair * Mesh_Manager::get_mesh_pair_pointer(moris_index aPairIndex)
        {
            MORIS_ASSERT( aPairIndex < (moris_index) mMeshPairs.size(),
                          "Mesh_Manager::get_mesh_pair: requested mesh pair does not exist.");
            return &mMeshPairs(aPairIndex);
        }

        //-------------------------------------------------------------------------

        void
        Mesh_Manager::get_mesh_pair(
                moris_index          aPairIndex,
                Interpolation_Mesh * & aInterpolationMesh,
                Integration_Mesh   * & aIntegrationMesh)
        {
            MORIS_ASSERT( aPairIndex < (moris_index) mMeshPairs.size(),
                    "Mesh_Manager::get_mesh_pair: requested mesh pair does not exist.");

            aInterpolationMesh = mMeshPairs(aPairIndex).get_interpolation_mesh();
            aIntegrationMesh = mMeshPairs(aPairIndex).get_integration_mesh();
        }

        //-------------------------------------------------------------------------

        Interpolation_Mesh*
        Mesh_Manager::get_interpolation_mesh(moris_index aMeshIndex)
        {
            MORIS_ASSERT( aMeshIndex < (moris_index) mMeshPairs.size(),
                    "Mesh_Manager::get_interpolation_mesh: requested mesh does not exist.");

            return mMeshPairs(aMeshIndex).get_interpolation_mesh();
        }

        //-------------------------------------------------------------------------

        Integration_Mesh*
        Mesh_Manager::get_integration_mesh(moris_index aMeshIndex)
        {
            MORIS_ASSERT( aMeshIndex < (moris_index) mMeshPairs.size(),
                    "Mesh_Manager::get_integration_mesh: requested mesh does not exist.");

            return mMeshPairs(aMeshIndex).get_integration_mesh();
        }

        //-------------------------------------------------------------------------
    }
}

