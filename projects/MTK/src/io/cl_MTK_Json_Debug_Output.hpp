//
// Created by frank on 12/4/23.
//

#ifndef MORIS_CL_MTK_JSON_DEBUG_OUTPUT_HPP
#define MORIS_CL_MTK_JSON_DEBUG_OUTPUT_HPP

#include <unordered_set>
#include "unordered_map"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_Json_Object.hpp"

namespace moris
{
    namespace mtk
    {
        class Json_Debug_Output
        {
          public:
            explicit Json_Debug_Output( Integration_Mesh const *aMesh )
                    : mMesh( aMesh ){};

            void write_to_json( std::string const &aFileName );

            void set_ig_vertex_displacements( std::unordered_map< moris::moris_index, Vector< moris::real > > const &aIGVertexDisplacements )
            {
                mIGVertexDisplacements = aIGVertexDisplacements;
            }

            void set_ip_vertex_displacements( std::unordered_map< moris::moris_index, Vector< moris::real > > const &aIPVertexDisplacements )
            {
                mIPVertexDisplacements = aIPVertexDisplacements;
            }

          private:
            Integration_Mesh const                                         *mMesh;
            std::string                                                     mFileName;
            std::unordered_set< Vertex const * >                            mAllIGVertices;
            std::unordered_set< Vertex const * >                            mAllIPVertices;
            std::unordered_set< Cell const * >                              mAllIGCells;
            std::unordered_set< Cell const * >                              mAllIPCells;
            std::unordered_map< moris::moris_index, Vector< moris::real > > mIGVertexDisplacements;
            std::unordered_map< moris::moris_index, Vector< moris::real > > mIPVertexDisplacements;

            Json serialize_side_sets( Vector< moris::mtk::Side_Set * > &aSideSets );
            Json serialize_side_clusters( Vector< Cluster const * > &aSideClusters );
            Json serialize_side_cluster( Cluster const *tCluster );
            Json serialize_double_side_sets( Vector< moris::mtk::Double_Side_Set * > &aDoubleSideSets );
            Json serialize_block_sets( Vector< moris::mtk::Block_Set * > &tBlockSets );
            Json serialize_cell_cluster( Cluster const *aCluster );
            Json serialize_cell_clusters( Vector< Cluster const * > &tClusters );
            Json serialize_all_vertices( const std::unordered_set< Vertex const * > &aVertices, const std::unordered_map< moris::moris_index, Vector< moris::real > > &aVertexDisplacements );
            Json serialize_all_cells( const std::unordered_set< Cell const * > &aCells );
            Json serialize_ig_cells( Vector< moris::mtk::Cell const * > &aCells );
            Json serialize_ip_cell( moris::mtk::Cell const *&aCell );
        };
    }    // namespace mtk
}    // namespace moris

#endif    // MORIS_CL_MTK_JSON_DEBUG_OUTPUT_HPP
