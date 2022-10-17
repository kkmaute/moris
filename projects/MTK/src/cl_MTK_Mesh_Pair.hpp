/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Pair.hpp
 *
 */

#ifndef MORIS_CL_MTK_MESH_PAIR_HPP
#define MORIS_CL_MTK_MESH_PAIR_HPP

#include "cl_MTK_Mesh_Manager.hpp"

namespace moris
{
    namespace mtk
    {
        class Interpolation_Mesh;
        class Integration_Mesh;

        class Mesh_Pair
        {

          private:
            Interpolation_Mesh* mInterpolationMesh;
            Integration_Mesh*   mIntegrationMesh;
            bool                mIsOwned;

          public:
            /**
             * Constructor
             *
             * @param aInterpolationMesh Interpolation mesh
             * @param aIntegrationMesh Integration mesh
             * @param aIsOwned If mesh pointers are owned by this mesh pair
             */
            Mesh_Pair(
                    Interpolation_Mesh* aInterpolationMesh,
                    Integration_Mesh*   aIntegrationMesh,
                    bool                aIsOwned = false );

            /**
             * Explicit copy constructor to ensure pointers are only owned once.
             *
             * @param aMeshPair Input mesh pair for copying
             */
            Mesh_Pair( const Mesh_Pair& aMeshPair );

            /**
             * Destructor (deletes mesh pointers if owned)
             */
            ~Mesh_Pair();

            /**
             * Gets the interpolation mesh in this mesh pair.
             *
             * @return Interpolation mesh
             */
            Interpolation_Mesh* get_interpolation_mesh() const;

            /**
             * Gets the integration mesh in this mesh pair.
             *
             * @return Integration mesh
             */
            Integration_Mesh* get_integration_mesh() const;

            /**
             * Copy constructor
             */
            const Mesh_Pair&
            operator=( const Mesh_Pair& aMeshPair )
            {
                mInterpolationMesh = aMeshPair.mInterpolationMesh;
                mIntegrationMesh   = aMeshPair.mIntegrationMesh;
                mIsOwned           = aMeshPair.mIsOwned;

                return *this;
            }

            /**
             * Friend functions in the mesh manager so it can properly transfer pointer ownership, without revealing
             * this to the outside.
             */
            friend uint Mesh_Manager::register_mesh_pair( Mesh_Pair& aMeshPair );
            friend uint Mesh_Manager::register_mesh_pair(
                    Interpolation_Mesh* aInterpolationMesh,
                    Integration_Mesh*   aIntegrationMesh,
                    bool                aIsOwned,
                    std::string const & aMeshPairName );
        };
    }    // namespace mtk
}    // namespace moris

#endif    // MORIS_CL_MTK_MESH_PAIR_HPP
