/*
 * cl_XTK_Enriched_Interpolation_Mesh.hpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_

#include "cl_Cell.hpp"

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_XTK_Interpolation_Vertex_Unzipped.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Vertex_Enrichment.hpp"
#include "cl_XTK_Enrichment.hpp"

using namespace moris;

namespace xtk
{

class Model;
class Enriched_Interpolation_Mesh: public mtk::Interpolation_Mesh
{
public:
    Enriched_Interpolation_Mesh(Model* aXTKModel);

    // MTK Mesh Core Functionality
    MeshType           get_mesh_type() const;
    moris::uint        get_spatial_dim() const;
    uint               get_num_entities( enum EntityRank aEntityRank ) const;
    Matrix< IndexMat > get_entity_connected_to_entity_loc_inds(moris_index aEntityIndex, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank) const;
    Matrix< IndexMat > get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const;
    moris_id           get_glb_entity_id_from_entity_loc_index(moris_index aEntityIndex,enum EntityRank aEntityRank) const;
    moris_index        get_loc_entity_ind_from_entity_glb_id(moris_id        aEntityId, enum EntityRank aEntityRank) const;
    Cell<mtk::Vertex const *> get_all_vertices() const;

    Matrix<IdMat>
    get_entity_connected_to_entity_glob_ids( moris_id     aEntityId,
                                             enum EntityRank aInputEntityRank,
                                             enum EntityRank aOutputEntityRank) const;


    /*
     * Convenient helper functions
     */
    Matrix<IdMat> convert_indices_to_ids(Matrix<IndexMat> const & aIndices,
                                         enum EntityRank          aEntityRank) const;

    Matrix<IndexMat> convert_ids_to_indices(Matrix<IdMat> const & aIds,
                                            enum EntityRank       aEntityRank) const;

    Matrix< DDRMat >
    get_node_coordinate( moris_index aNodeIndex ) const;

    mtk::Vertex &       get_mtk_vertex( moris_index aVertexIndex );
    mtk::Vertex const & get_mtk_vertex( moris_index aVertexIndex ) const;
    mtk::Cell &         get_writable_mtk_cell( moris_index aElementIndex );

    Matrix< IdMat >     get_communication_table() const;

    /*
     * Number of interpolation cells
     */
    uint
    get_num_elements();

    uint
    get_num_verts_per_interp_cell();

    Interpolation_Vertex_Unzipped*
    get_unzipped_vertex_pointer(moris_index aVertexIndex);


    Cell<Interpolation_Cell_Unzipped> const &
    get_enriched_interpolation_cells() const;

    // Print functions
    void print() const;
    void print_enriched_cells() const;
    void print_vertex_maps() const;
    void print_enriched_cell_maps() const;


    // friend class
    friend class Enrichment;
protected:
    // Model pointer
    Model* mXTKModel;

    // enriched interpolation vertices
    moris::uint                         mNumVerts;
    Cell<Interpolation_Vertex_Unzipped> mEnrichedInterpVerts;

    // enriched interpolation cells
    moris::uint                         mNumVertsPerInterpCell;
    Cell<Interpolation_Cell_Unzipped>   mEnrichedInterpCells;

    // for each outer cell (base interpolation vertex), pointer to the various interpolation vertex enrichment
    Cell<Cell<moris_index>> mBaseInterpVertToVertEnrichmentIndex;
    Cell<Vertex_Enrichment> mInterpVertEnrichment;
    Cell<moris_index>       mVertexEnrichmentParentVertexIndex;

    // cell mapping
    Cell<Matrix< IdMat >>                          mLocalToGlobalMaps;
    Cell<std::unordered_map<moris_id,moris_index>> mGlobaltoLobalMaps;

    // functions used by enrichment for construction of the mesh
    /*
     * Add a vertex enrichment to the member data. returns the index of the vertex enrichment.
     * I do not return a pointer here because addresses may change while constructing the data
     */
    moris_index
    add_vertex_enrichment( mtk::Vertex *       aBaseInterpVertex,
                           Vertex_Enrichment & aVertexEnrichment,
                           bool              & aNewVertex);

    /*
     * Get the pointer to the vertex enrichment provided the vertex enrichment index.
     */
    Vertex_Enrichment*
    get_vertex_enrichment(moris_index aVertexEnrichmentIndex);

    /*
     * Returns the vertex index corresponding to the vertex enrichment
     */
    moris_index
    get_vertex_related_to_vertex_enrichment(moris_index aVertexEnrichmentIndex) const;

    void
    finalize_setup();

    // map setup
    void setup_local_to_global_maps();
    void setup_vertex_maps();
    void setup_cell_maps();


};


}


#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_ */
