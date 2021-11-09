/**
 * cl_XTK_Elevate_Order_Interface.hpp  
 * 
 *  Created on: Nov  03, 2021 
 *      Author: Nils Wunsch
 */
#ifndef SRC_cl_XTK_Elevate_Order_Interface
#define SRC_cl_XTK_Elevate_Order_Interface

#include "cl_XTK_Decomposition_Algorithm.hpp"

namespace moris
{
namespace mtk
{
    class Mesh;
}
}// namespace moris

namespace xtk
{

class Integration_Mesh_Generator;
class Integration_Mesh_Generation_Data;
class Decomposition_Data;
class Cut_Integration_Mesh;
class IG_Cell_Group;

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

// pure virtual base class for 2D and 3D Regular Subdivision Data
class Elevate_Order_Template
{
  public:
    Elevate_Order_Template() {}

    virtual ~Elevate_Order_Template() {}

    /**
     * @brief Get number of nodes to be added per element
     * 
     * @return moris_index number of nodes to be added per element
     */
    virtual moris_index
    get_num_new_nodes() const = 0;

    /**
     * @brief Get total number of IG vertices on new elements for Template
     * 
     * @return moris_index total number of IG vertices on new elements
     */
    virtual moris_index
    get_total_ig_verts() const = 0;

    /**
     * @brief gives back whether new vertices will be created on the requested entity rank
     * 
     * @param aEntityRank Entity Rank for which new vertices should be created
     * @return whether new vertices will be created or not
     */
    virtual bool
    has_new_vertices_on_entity( enum EntityRank aEntityRank ) const = 0;

    /**
     * @brief gives back the number of new vertices that will be created on the requested entity rank
     * 
     * @param aEntityRank Entity Rank for which new vertices should be created
     * @return number of new vertices on that entity rank
     */
    virtual uint
    num_new_vertices_per_entity( enum EntityRank aEntityRank ) const = 0;

    /**
     * @brief Get parametric coords for new vertices wrt entity they are created on
     * 
     * @param aEntityRank Entity Rank for which new vertices should be created
     * @return Cell< Matrix< DDRMat > > list of parametric coords for vertices to be created per entity
     */
    virtual Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coords_wrt_entity( enum EntityRank aEntityRank ) const = 0;

    /**
     * @brief Get the new ig cell's topology
     * 
     * @return enum CellTopology new ig cell topology
     */
    virtual enum CellTopology
    get_ig_cell_topology() const = 0;

    /**
     * @brief Get number of spatial dimensions current template assumes
     * 
     * @return uint number of spatial dimensions
     */
    virtual uint
    get_num_spatial_dims() const = 0;
};

// -------------------------------------------------------------------------

class TRI3_to_TRI6 : public Elevate_Order_Template
{
  public:
    TRI3_to_TRI6() {}

    moris_index
    get_num_new_nodes() const
    {
        return 3;
    }

    moris_index
    get_total_ig_verts() const
    {
        return 6;
    }

    bool
    has_new_vertices_on_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::has_new_nodes_on_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return false;
            break;
        }
        case EntityRank::EDGE:
        {
            return true;
            break;
        }
        case EntityRank::FACE:
        {
            return false;
            break;
        }
        case EntityRank::ELEMENT:
        {
            return false;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::has_new_nodes_on_entity() - Unknown Entity Rank." );
            return false;
            break;
        }
        }
    }

    uint
    num_new_vertices_per_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::num_new_vertices_per_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return 0;
            break;
        }
        case EntityRank::EDGE:
        {
            return 1;
            break;
        }
        case EntityRank::FACE:
        {
            return 0;
            break;
        }
        case EntityRank::ELEMENT:
        {
            return 0;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::num_new_vertices_per_entity() - Unknown Entity Rank." );
            return 0;
            break;
        }
        }
    }

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coords_wrt_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        case EntityRank::EDGE:
        {
            return { { { 0.0 } } };
            break;
        }
        case EntityRank::FACE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on faces. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        case EntityRank::ELEMENT:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - No new vertices inside cell. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TRI3_to_TRI6::get_new_vertex_parametric_coords_wrt_entity() - Unknown Entity Rank." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        }
    }

    enum CellTopology
    get_ig_cell_topology() const
    {
        return CellTopology::TRI6;
    }

    uint
    get_num_spatial_dims() const
    {
        return 2;
    }
};

// -------------------------------------------------------------------------

class TET4_to_TET10 : public Elevate_Order_Template
{
  public:
    TET4_to_TET10() {}

    moris_index
    get_num_new_nodes() const
    {
        return 6;
    }

    moris_index
    get_total_ig_verts() const
    {
        return 10;
    }

    bool
    has_new_vertices_on_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TET4_to_TET10::has_new_nodes_on_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return false;
            break;
        }
        case EntityRank::EDGE:
        {
            return true;
            break;
        }
        case EntityRank::FACE:
        {
            return false;
            break;
        }
        case EntityRank::ELEMENT:
        {
            return true;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TET4_to_TET10::has_new_nodes_on_entity() - Unknown Entity Rank." );
            return false;
            break;
        }
        }
    }

    uint
    num_new_vertices_per_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TET4_to_TET10::num_new_vertices_per_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return 0;
            break;
        }
        case EntityRank::EDGE:
        {
            return 1;
            break;
        }
        case EntityRank::FACE:
        {
            return 0;
            break;
        }
        case EntityRank::ELEMENT:
        {
            return 1;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TET4_to_TET10::num_new_vertices_per_entity() - Unknown Entity Rank." );
            return 0;
            break;
        }
        }
    }

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coords_wrt_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        case EntityRank::EDGE:
        {
            return { { { 0.0 } } };
            break;
        }
        case EntityRank::FACE:
        {
            MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on faces. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        case EntityRank::ELEMENT:
        {
            MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices inside cell. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TET4_to_TET10::get_new_vertex_parametric_coords_wrt_entity() - Unknown Entity Rank." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        }
    }

    enum CellTopology
    get_ig_cell_topology() const
    {
        return CellTopology::TET10;
    }

    uint
    get_num_spatial_dims() const
    {
        return 3;
    }
};

// -------------------------------------------------------------------------

class TRI3_to_TRI10 : public Elevate_Order_Template
{
  public:
    TRI3_to_TRI10() {}

    moris_index
    get_num_new_nodes() const
    {
        return 7;
    }

    moris_index
    get_total_ig_verts() const
    {
        return 10;
    }

    bool
    has_new_vertices_on_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::has_new_nodes_on_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return false;
            break;
        }
        case EntityRank::EDGE:
        {
            return true;
            break;
        }
        case EntityRank::FACE:
        {
            return false;
            break;
        }
        case EntityRank::ELEMENT:
        {
            return false;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::has_new_nodes_on_entity() - Unknown Entity Rank." );
            return false;
            break;
        }
        }
    }

    uint
    num_new_vertices_per_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::num_new_vertices_per_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return 0;
            break;
        }
        case EntityRank::EDGE:
        {
            return 2;
            break;
        }
        case EntityRank::FACE:
        {
            return 0;
            break;
        }
        case EntityRank::ELEMENT:
        {
            return 1;
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::num_new_vertices_per_entity() - Unknown Entity Rank." );
            return 0;
            break;
        }
        }
    }

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coords_wrt_entity( enum EntityRank aEntityRank ) const
    {
        switch ( aEntityRank )
        {
        case EntityRank::NODE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on vertices themselves. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        case EntityRank::EDGE:
        {
            return {
                { { -1.0 / 3.0 } },
                { { 1.0 / 3.0 } }
            };
            break;
        }
        case EntityRank::FACE:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::get_new_vertex_parametric_coords_wrt_entity() - No new vertices on faces. This shouldn't be requested." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        case EntityRank::ELEMENT:
        {
            real tOneThird = 1.0 / 3.0;
            return { { { tOneThird, tOneThird, tOneThird } } };
            break;
        }
        default:
        {
            MORIS_ERROR( false, "TRI3_to_TRI10::get_new_vertex_parametric_coords_wrt_entity() - Unknown Entity Rank." );
            return Cell< Matrix< DDRMat > >( 0 );
            break;
        }
        }
    }

    enum CellTopology
    get_ig_cell_topology() const
    {
        return CellTopology::TRI10;
    }

    uint
    get_num_spatial_dims() const
    {
        return 2;
    }
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

class Elevate_Order_Interface : public Decomposition_Algorithm
{

    // -------------------------------------------------------------------------

  private:
    std::shared_ptr< Elevate_Order_Template > mElevateOrderTemplate;

    Integration_Mesh_Generation_Data* mMeshGenerationData;
    Decomposition_Data*               mDecompositionData;
    Cut_Integration_Mesh*             mCutIntegrationMesh;
    moris::mtk::Mesh*                 mBackgroundMesh;
    Integration_Mesh_Generator*       mGenerator;
    moris::uint                       mNumTotalCells = 0;

    // -------------------------------------------------------------------------

  public:
    Elevate_Order_Interface( ParameterList& aParameterList, enum Subdivision_Method aSubdivisionMethod );

    ~Elevate_Order_Interface() {}

    moris_index get_signature() const;

    enum Decomposition_Algorithm_Type get_algorithm_type() const;

    bool has_geometric_independent_vertices() const;

    void perform(
        Integration_Mesh_Generation_Data* aMeshGenerationData,
        Decomposition_Data*               aDecompositionData,
        Cut_Integration_Mesh*             aCutIntegrationMesh,
        moris::mtk::Mesh*                 aBackgroundMesh,
        Integration_Mesh_Generator*       aMeshGenerator );

    void
    perform_impl_vertex_requests(
        Integration_Mesh_Generation_Data* aMeshGenerationData,
        Decomposition_Data*               aDecompositionData,
        Cut_Integration_Mesh*             aCutIntegrationMesh,
        moris::mtk::Mesh*                 aBackgroundMesh,
        Integration_Mesh_Generator*       aMeshGenerator );

    void
    perform_impl_generate_mesh(
        Integration_Mesh_Generation_Data* aMeshGenerationData,
        Decomposition_Data*               aDecompositionData,
        Cut_Integration_Mesh*             aCutIntegrationMesh,
        moris::mtk::Mesh*                 aBackgroundMesh,
        Integration_Mesh_Generator*       aMeshGenerator );

    // -------------------------------------------------------------------------

  private:
    bool
    make_vertex_requests(
        std::shared_ptr< Edge_Based_Connectivity >         aEdgeConnectivity,
        std::shared_ptr< Edge_Based_Ancestry >             aIgEdgeAncestry,
        moris::Cell< moris::mtk::Cell* >*                  aBackgroundCellForEdge,
        moris::Cell< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups,
        moris::Cell< moris::mtk::Cell* >*                  aIgCells );

    bool
    associate_new_vertices_with_cell_groups(
        std::shared_ptr< Edge_Based_Connectivity >         aEdgeConnectivity,
        std::shared_ptr< Edge_Based_Ancestry >             aIgEdgeAncestry,
        moris::Cell< moris::mtk::Cell* >*                  aBackgroundCellForEdge,
        moris::Cell< std::shared_ptr< IG_Vertex_Group > >* aVertexGroups,
        moris::Cell< moris::mtk::Cell* >*                  aIgCells );

    void
    create_higher_order_integration_cells(
        std::shared_ptr< Edge_Based_Connectivity > aEdgeConnectivity,
        std::shared_ptr< Edge_Based_Ancestry >     aIgEdgeAncestry,
        moris::Cell< moris::mtk::Cell* >*          aIgCells );

    // -------------------------------------------------------------------------

    /**
     * @brief Creates unique id for an edge based on the IDs of its two end vertices
     * using the Cantor pairing function
     * 
     * @param aEdgeVertices list of mtk::vertices on edge
     * @return moris_index unique id for edge
     */
    moris_index
    hash_edge( moris::Cell< moris::mtk::Vertex* > const& aEdgeVertices );

    /**
     * @brief Creates unique id for a face based on the IDs of its three end corner vertices
     * by recursively applying the Cantor pairing function
     * 
     * @param aFaceVertices list of mtk::vertices at face corners
     * @return moris_index unique id for face
     */
    moris_index
    hash_face( moris::Cell< moris::mtk::Vertex* > const& aFaceVertices );

    /**
     * @brief swaps the values of two inidices
     * 
     * @param aInd1 first index, value to be swaped with ...
     * @param aInd2 the second index
     */
    void
    swap_indices( moris_index& aInd1, moris_index& aInd2 );

    // -------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_edge_vertex_global_coordinates( moris::Cell< moris::mtk::Vertex* > const& aEdgeVertices, Matrix< DDRMat > aEdgeCoord );

    Matrix< DDRMat >
    compute_tri_vertex_global_coordinates( moris::Cell< moris::mtk::Vertex* > const& aTriVertices, Matrix< DDRMat > aTriCoords );

    Matrix< DDRMat >
    compute_tet_vertex_global_coordinates( moris::Cell< moris::mtk::Vertex* > const& aTetVertices, Matrix< DDRMat > aTetCoords );
};


}// namespace xtk
#endif /* cl_XTK_Elevate_Order_Interface.hpp */