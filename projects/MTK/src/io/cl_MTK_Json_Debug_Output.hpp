//
// Created by frank on 12/4/23.
//

#ifndef MORIS_CL_MTK_JSON_DEBUG_OUTPUT_HPP
#define MORIS_CL_MTK_JSON_DEBUG_OUTPUT_HPP

#include <unordered_set>
#include "json.hpp"
#include "unordered_map"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"

using json = nlohmann::json;

namespace moris
{
    namespace mtk
    {
        class Json_Debug_Output
        {
          public:
            Json_Debug_Output( Integration_Mesh_DataBase_IG const *aMesh )
                    : mMesh( aMesh ){};

            void write_to_json( std::string const &aFileName );

          private:
            Integration_Mesh_DataBase_IG const  *mMesh;
            std::string                          mFileName;
            std::unordered_set< Vertex const * > mAllIGVertices;
            std::unordered_set< Vertex const * > mAllIPVertices;
            std::unordered_set< Cell const * >   mAllIGCells;
            std::unordered_set< Cell const * >   mAllIPCells;

            json          serialize_side_sets( moris::Cell< moris::mtk::Side_Set          *> &aSideSets );
            json::array_t serialize_side_clusters( moris::Cell< Cluster const * > &aSideClusters );
            json          serialize_side_cluster( Cluster const *tCluster );
            json          serialize_double_side_sets( moris::Cell< moris::mtk::Double_Side_Set          *> &aDoubleSideSets );
            json          serialize_block_sets( moris::Cell< moris::mtk::Block_Set          *> &tBlockSets );
            json          serialize_cell_cluster( Cluster const *aCluster );
            json::array_t serialize_cell_clusters( moris::Cell< Cluster const * > &tClusters );
            json          serialize_all_vertices( const std::unordered_set< Vertex const          *> &aVertices );

            json::array_t serialize_all_cells( const std::unordered_set< Cell const * >& aCells );
            json::array_t serialize_cells( moris::Cell< moris::mtk::Cell const * > &aCells, bool aIsIPCell = false );
        };
    }    // namespace mtk
}    // namespace moris

#endif    // MORIS_CL_MTK_JSON_DEBUG_OUTPUT_HPP
