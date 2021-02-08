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
            Interpolation_Mesh* mInterpolationMesh; //! Interpolation mesh
            Integration_Mesh* mIntegrationMesh; //! Integration mesh
            bool mIsOwned = false; //! If the mesh pointers are owned by the mesh pair

            /**
             * Destructor (deletes mesh pointers if owned)
             */
            ~Mesh_Pair();
        };
    }
}

#endif //MORIS_ST_MTK_MESH_PAIR_HPP
