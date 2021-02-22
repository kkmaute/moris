#ifndef MORIS_ST_MTK_MESH_PAIR_HPP
#define MORIS_ST_MTK_MESH_PAIR_HPP

namespace moris
{
    namespace mtk
    {
        class Interpolation_Mesh;
        class Integration_Mesh;

        struct Mesh_Pair
        {
            Interpolation_Mesh* mInterpolationMesh = nullptr; //! Interpolation mesh
            Integration_Mesh* mIntegrationMesh = nullptr; //! Integration mesh
            bool mIsOwned = false; //! If the mesh pointers are owned by the mesh pair

            /**
             * Default constructor
             */
            Mesh_Pair() = default;

            /**
             * Explicit copy constructor to ensure pointers are only owned once.
             *
             * @param aMeshPair Input mesh pair for copying
             */
            Mesh_Pair(const Mesh_Pair& aMeshPair);

            /**
             * Destructor (deletes mesh pointers if owned)
             */
            ~Mesh_Pair();
        };
    }
}

#endif //MORIS_ST_MTK_MESH_PAIR_HPP
