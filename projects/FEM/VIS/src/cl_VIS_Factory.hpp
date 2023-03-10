/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_VIS_FACTORY_HPP_
#define SRC_FEM_CL_VIS_FACTORY_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Vertex_Visualization.hpp"
#include "cl_VIS_Cell_Visualization.hpp"
#include "cl_VIS_Cell_Cluster_Visualization.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Output_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        class Set;
        class Cluster;
        class Cell;
        class Vertex;
        class Mesh_Manager;
    }

    namespace vis
    {
        struct Output_Data;

        class VIS_Factory
        {
        private:
            moris::Cell< moris::mtk::Set * > mListOfBlocks;

            moris::Cell< moris::Cell< const mtk::Cluster * > >        mClustersOnBlock;   //FIXME delete can be used temporary
            moris::Cell< moris::Cell< mtk::Cell * > >   mCellsOnSet;
            moris::Cell< moris::Cell< mtk::Vertex * > > mVerticesOnSet;

            moris::Cell< Matrix< DDSMat > >             mVertexMapOnSet;
            moris::Cell< Matrix< DDSMat > >             mCellMapOnSet;

            uint mNumRequestedSets;
            bool mOnlyPrimaryCells = false;

            std::shared_ptr< mtk::Mesh_Manager > mMesh = nullptr;
            const uint         mMeshPairIndex;

            moris::Cell< std::string > tRequestedSetNames;

        public:
            VIS_Factory(
                    std::shared_ptr< mtk::Mesh_Manager > aMesh,
                    const uint         aMeshPairIndex ) : mMesh( aMesh ),
                                                           mMeshPairIndex( aMeshPairIndex )
            {};

//-----------------------------------------------------------------------------------------------------------

            ~VIS_Factory(){};

//-----------------------------------------------------------------------------------------------------------

            mtk::Mesh * create_visualization_mesh( moris::vis::Output_Data & aOutputData );

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_vertices();

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_cells();

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_clusters( moris::vis::Output_Data & aOutputData );

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_blocks();

//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_FACTORY_HPP_ */

