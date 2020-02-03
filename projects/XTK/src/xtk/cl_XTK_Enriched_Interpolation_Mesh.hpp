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
#include "cl_XTK_Ghost_Stabilization.hpp"
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
    ~Enriched_Interpolation_Mesh();
    //------------------------------------------------------------------------------
    // MTK Mesh Core Impl (see base class mtk::Mesh for function descriptions)
    //------------------------------------------------------------------------------
    MeshType                  get_mesh_type() const;
    moris::uint               get_spatial_dim() const;
    uint                      get_num_entities( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
    uint                      get_num_coeffs(const uint aBSplineMeshIndex) const;
    Matrix< IndexMat >        get_entity_connected_to_entity_loc_inds(moris_index aEntityIndex, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank, const moris_index aIndex = 0) const;
    Matrix< IndexMat >        get_elements_connected_to_element_and_face_ind_loc_inds(moris_index aElementIndex) const;
    moris_id                  get_glb_entity_id_from_entity_loc_index(moris_index aEntityIndex,enum EntityRank aEntityRank,const moris_index aIndex = 0) const;
    moris_index               get_loc_entity_ind_from_entity_glb_id( moris_id aEntityId, enum EntityRank aEntityRank,const moris_index aIndex = 0) const;
    Cell<mtk::Vertex const *> get_all_vertices() const;
    Matrix<IdMat>             get_entity_connected_to_entity_glob_ids( moris_id aEntityId, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank, const moris_index aIndex = 0) const;
    Matrix< DDRMat >          get_node_coordinate( moris_index aNodeIndex ) const;
    mtk::Vertex &             get_mtk_vertex( moris_index aVertexIndex );
    mtk::Vertex const &       get_mtk_vertex( moris_index aVertexIndex ) const;
    mtk::Cell   const &       get_mtk_cell( moris_index aElementIndex ) const;
    mtk::Cell &               get_writable_mtk_cell( moris_index aElementIndex );
    Matrix< IdMat >           get_communication_table() const;
    uint                      get_num_elements();
    moris_id                  get_max_entity_id( enum EntityRank aEntityRank,const moris_index aIndex =0 ) const;
    void                      get_adof_map( const uint aBSplineIndex, map< moris_id, moris_index > & aAdofMap ) const;

    //------------------------------------------------------------------------------
    // end mesh core functions
    //------------------------------------------------------------------------------

    //------------------------------------------------------------------------------
    // MTK Interpolation Mesh Impl
    // see base class mtk::Interpolation_Mesh for function descriptions
    //------------------------------------------------------------------------------
    // PLACEHOLDER
    //------------------------------------------------------------------------------
    // end interpolation mesh functions
    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // Accessor functions of XTK specific data structures
    //------------------------------------------------------------------------------
    /*!
     * Get the enriched interpolation coefficients associated with a background coefficient
     */
    Matrix<IndexMat> const &
    get_enriched_coefficients_at_background_coefficient(moris_index aBackgroundCoeffIndex) const;
    //------------------------------------------------------------------------------
    /*!
     * get the enriched coefficients at all the background mesh coefficients
     */
    Cell<Matrix<IndexMat>> const &
    get_enriched_coefficients_to_background_coefficients() const;
    //------------------------------------------------------------------------------
    /*!
     * get the local enriched coefficient to global map
     */
    Matrix<IndexMat> const &
    get_enriched_coefficient_local_to_global_map() const;
    //------------------------------------------------------------------------------
    /*!
     * Return the vector of background coefficient local to global
     */
    Matrix<IndexMat>
    get_background_coefficient_local_to_global_map() const;
    //------------------------------------------------------------------------------
    uint
    get_num_background_coefficients() const;

    /*!
     * Returns the number of vertices per interpolation cell
     */
    uint
    get_num_verts_per_interp_cell();
    //------------------------------------------------------------------------------
    /*
     * Returns the interpolation vertex unzipped for provided vertex index
     */
    Interpolation_Vertex_Unzipped*
    get_unzipped_vertex_pointer(moris_index aVertexIndex);
    //------------------------------------------------------------------------------
    /*!
     * Return the enriched interpolation cells
     */
    Cell<Interpolation_Cell_Unzipped*> const &
    get_enriched_interpolation_cells() const;

protected:
    Cell<Interpolation_Cell_Unzipped*> &
    get_enriched_interpolation_cells() ;

    /*!
     * Return the owned interpolation cells, shared cells sorted by owning proc, owning procs
     */
    void
    get_owned_and_not_owned_enriched_interpolation_cells(Cell<Interpolation_Cell_Unzipped *>      & aOwnedInterpCells,
                                                        Cell<Cell<Interpolation_Cell_Unzipped *>> & aNotOwnedInterpCells,
                                                        Cell<uint>                                & aProcRanks);
public:
    //------------------------------------------------------------------------------
    /*
     * Provided pointers to base interpolation cells, return all the enriched interpolation cells attached to these
     * cells
     */
    moris::Cell<Interpolation_Cell_Unzipped const*>
    get_enriched_cells_from_base_cells(moris::Cell<moris::mtk::Cell const*> const & aBaseCells) const;
    //------------------------------------------------------------------------------
    /*!
     *  Single cell version of above
     */
    moris::Cell<Interpolation_Cell_Unzipped const*>
    get_enriched_cells_from_base_cell(moris::mtk::Cell const * aBaseCells) const;
    //------------------------------------------------------------------------------
    Interpolation_Vertex_Unzipped const &
    get_xtk_interp_vertex(moris::uint aVertexIndex) const;
    //------------------------------------------------------------------------------
    /*
     * Convert a entity indices to entity ids
     */
    Matrix<IdMat> convert_indices_to_ids(Matrix<IndexMat> const & aIndices,
                                         enum EntityRank          aEntityRank) const;
    //------------------------------------------------------------------------------
    /*
     * Convert a entity ids to entity indices
     */
    Matrix<IndexMat> convert_ids_to_indices(Matrix<IdMat> const & aIds,
                                            enum EntityRank       aEntityRank) const;
    //------------------------------------------------------------------------------
    /*!
     * convert enriched
     */
    //------------------------------------------------------------------------------
    void
    convert_enriched_basis_indices_to_ids(Matrix<IndexMat> const & aEnrichedIndices,
                                          Matrix<IdMat>          & aEnrichedIds) const;

    void
    convert_enriched_basis_indices_to_ids(Cell<Matrix<IndexMat>> const & aEnrichedIndices,
                                          Cell<Matrix<IdMat>>          & aEnrichedIds) const;

    // Print functions
    void print() const;
    void print_enriched_cells() const;
    void print_vertex_maps() const;
    void print_enriched_cell_maps() const;
    void print_basis_to_enriched_basis() const;
    void print_vertex_interpolation() const;


    // friend class
    friend class Enrichment;
    friend class Enriched_Integration_Mesh;
    friend class Ghost_Stabilization;
protected:
    // Model pointer
    Model* mXTKModel;

    // basis rank
    enum EntityRank mBasisRank;

    // mesh index
    moris_index mMeshIndex;

    // enriched interpolation vertices
    moris::uint                         mNumVerts;
    Cell<Interpolation_Vertex_Unzipped> mEnrichedInterpVerts;

    // enriched interpolation cells
    moris::uint                         mNumVertsPerInterpCell;
    Cell<Interpolation_Cell_Unzipped*>  mEnrichedInterpCells;

    // for each outer cell (base interpolation vertex), indices of enriched vertices
    Cell<Cell<moris_index>> mBaseInterpVertToVertEnrichmentIndex;

    // vertex enrichments
    Cell<Vertex_Enrichment *> mInterpVertEnrichment;

    // vertex enrichment to parent vertex index (these are enriched interpolation vertex indices)
    Cell<moris_index>       mVertexEnrichmentParentVertexIndex;

    // basis coefficient to enriched basis coefficient
    moris::Cell<moris::Matrix<moris::IndexMat>> mCoeffToEnrichCoeffs;

    // local to global enriched basis vector
    moris::Matrix<moris::IdMat> mEnrichCoeffLocToGlob;

    // Entity maps
    Cell<Matrix< IdMat >>                          mLocalToGlobalMaps;
    Cell<std::unordered_map<moris_id,moris_index>> mGlobaltoLobalMaps;

    // base interpolation cells to their enriched interpolation cells
    moris::Cell<moris::Cell<Interpolation_Cell_Unzipped*>> mBaseCelltoEnrichedCell;

    // a connecitivty pointer that all the enriched interpolation cells use
    moris::mtk::Cell_Info* mCellInfo;

    // functions used by enrichment for construction of the mesh
    /*
     * Add a vertex enrichment to the member data. returns the index of the vertex enrichment.
     * I do not return a pointer here because addresses may change while constructing the data
     */
    moris_index
    add_vertex_enrichment( mtk::Vertex *       aBaseInterpVertex,
                           Vertex_Enrichment & aVertexEnrichment,
                           bool              & aNewVertex);
    //------------------------------------------------------------------------------
    /*
     * Get the pointer to the vertex enrichment provided the vertex enrichment index.
     */
    Vertex_Enrichment*
    get_vertex_enrichment(moris_index aVertexEnrichmentIndex);
    //------------------------------------------------------------------------------
    /*
     * Returns the vertex index corresponding to the vertex enrichment
     */
    moris_index
    get_vertex_related_to_vertex_enrichment(moris_index aVertexEnrichmentIndex) const;

    void
    finalize_setup();

    void
    assign_vertex_interpolation_ids();

    // map setup
    void setup_local_to_global_maps();
    void setup_vertex_maps();
    void setup_cell_maps();

    //------------------------------------------------------------------------------
    // Parallel functions
    //------------------------------------------------------------------------------
    moris_id allocate_entity_ids(moris::size_t  aNumReqs, enum EntityRank aEntityRank);


};


}


#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_ */
