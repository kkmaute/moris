/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_Manager.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Map.hpp"

namespace moris
{
    namespace hmr
    {
        class HMR;
    }

    namespace mtk
    {
        class Mesh_Pair;
        class Interpolation_Mesh;
        class Integration_Mesh;
        class Field;

        class Mesh_Manager : public std::enable_shared_from_this< Mesh_Manager >
        {
        private:
            //! list of mesh pairs
            moris::Cell<Mesh_Pair>                     mMeshPairs;

            //! mesh pair to index map
            moris::map< std::string, moris_index >     mMeshPairNameToIndexMap;

            //! pointer to HMR performer
            std::weak_ptr<hmr::HMR>                    mHMRPerformer ;

            //! list of registered fields
            moris::Cell< std::weak_ptr< mtk::Field > > mFields;

            //! mesh pair to index map
            moris::map< std::string, moris_index >     mFieldLabelToIndexMap;

            moris::Cell< moris::Cell< moris_index > >  mMeshPairToFieldIndexMap;

        public:

            Mesh_Manager();

            //--------------------------------------------------------------------

            ~Mesh_Manager();

            //--------------------------------------------------------------------

            void set_performer( std::shared_ptr< hmr::HMR > aHMRPerformer )
            {
                mHMRPerformer = aHMRPerformer;
            };

            //--------------------------------------------------------------------

            uint register_mesh_pair(
                    Interpolation_Mesh* aInterpolationMesh,
                    Integration_Mesh*   aIntegrationMesh,
                    bool                aIsOwned = false,
                    std::string const & aMeshPairName = "" );

            //--------------------------------------------------------------------

            /**
             * Register a mesh pair with the mesh manager.
             *
             * @param aMeshPair Mesh pair
             * @return Mesh pair index
             */
            uint register_mesh_pair(Mesh_Pair& aMeshPair);

            //--------------------------------------------------------------------

            /**
             * Gets a mesh pair.
             *
             * @param aPairIndex Mesh pair index
             * @return Mesh pair
             */
            const Mesh_Pair & get_mesh_pair(moris_index aPairIndex);

            //--------------------------------------------------------------------

            /**
             * Gets a mesh pair.
             *
             * @param aPairIndex Mesh pair name
             * @return Mesh pair
             */
            const Mesh_Pair & get_mesh_pair( const std::string & aMeshPairName );

            //--------------------------------------------------------------------

            /**
             * Gremoves mesh pair with name
             *
             * @param aPairIndex Mesh pair name
             */
            void remove_mesh_pair( const std::string & aMeshPairName );

            //-------------------------------------------------------------------------

            void update_mesh_pairs( const moris::Cell< std::string > & aMeshPairNames );

            //--------------------------------------------------------------------

            void
            get_mesh_pair(
                    moris_index            aPairIndex,
                    Interpolation_Mesh * & aInterpolationMesh,
                    Integration_Mesh   * & aIntegrationMesh);

            //--------------------------------------------------------------------

            Interpolation_Mesh*
            get_interpolation_mesh(moris_index aMeshIndex);

            //--------------------------------------------------------------------

            Integration_Mesh*
            get_integration_mesh(moris_index aMeshIndex);

            //--------------------------------------------------------------------

            std::shared_ptr< Mesh_Manager > get_pointer()
            {
                return shared_from_this();
            }

            //--------------------------------------------------------------------

            void register_field( std::shared_ptr< mtk::Field > aField );

            //--------------------------------------------------------------------

            //--------------------------------------------------------------------

            //--------------------------------------------------------------------
        };
    }
}
#endif /* PROJECTS_MTK_SRC_CL_MTK_MESH_MANAGER_HPP_ */

