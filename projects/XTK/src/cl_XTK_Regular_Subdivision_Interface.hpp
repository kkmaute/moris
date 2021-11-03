#ifndef MORIS_CL_XTK_REGULAR_SUBDIVISION_INTERFACE_HPP_
#define MORIS_CL_XTK_REGULAR_SUBDIVISION_INTERFACE_HPP_

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
struct Regular_Subdivision_Interface_Data
{
  public:
    Matrix< IndexMat >       mNewNodesOnEdges;
    Matrix< IndexMat >       mNewNodesOnFaces;
    Matrix< IndexMat >       mNewNodesOnCells;
    Matrix< IndexMat >       mNewNodesOnEdgesOrd;
    Matrix< IndexMat >       mNewNodesOnFacesOrd;
    Matrix< IndexMat >       mNewNodesOnCellsOrd;
    Vertex_Ancestry          mVertexAncestry;
    Cell< Matrix< DDRMat > > mNewNodeXi;
    Matrix< DDRMat >         mNXi;
};

class Integration_Mesh_Generator;
class Integration_Mesh_Generation_Data;
class Decomposition_Data;
class Cut_Integration_Mesh;

// -------------------------------------------------------------------------

// pure virtual base class for 2D and 3D Regular Subdivision Data
class Regular_Subdivision_Template
{
  public:

    Regular_Subdivision_Template() {}

    virtual
    ~Regular_Subdivision_Template() {}

    virtual
    moris_index
    get_num_new_nodes() const = 0;

    virtual
    moris_index
    get_total_ig_verts() const = 0;

    virtual
    moris_index
    get_num_ig_cells() const = 0;

    virtual
    moris_index
    get_num_verts_per_cell() const = 0;

    virtual
    Matrix< IndexMat >
    get_new_node_on_parent_edge() const = 0;

    virtual
    Matrix< IndexMat >
    get_new_node_on_parent_edge_edge_ordinal() const = 0;

    virtual
    Matrix< IndexMat >
    get_new_node_on_parent_face() const = 0;

    virtual
    Matrix< IndexMat >
    get_new_node_on_parent_face_face_ordinal() const = 0;

    virtual
    Matrix< IndexMat >
    get_new_node_in_parent_cell() const = 0;

    virtual
    Matrix< IndexMat >
    get_new_node_on_parent_cell_cell_ordinal() const = 0;

    virtual
    Vertex_Ancestry
    get_vertex_ancestry() const = 0;

    virtual
    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coordinates_wrt_parent() const = 0;

    virtual
    Cell< Cell< moris_index > >
    get_ig_cell_to_vertex_connectivity() const = 0;

    virtual
    moris_index
    get_parametric_dimension() const = 0;

    virtual
    moris_index
    get_signature() const = 0;

    virtual
    enum CellTopology
    get_ig_cell_topology() const = 0;
};

// -------------------------------------------------------------------------

// 3 ----- 2
// | \   / |
// |   4   |
// | /   \ |
// 0 ----- 1

class Regular_Subdivision_4_TRIS : public Regular_Subdivision_Template
{
  public:
    Regular_Subdivision_4_TRIS() {}

    moris_index
    get_num_new_nodes() const
    {
        return 1;
    }

    moris_index
    get_total_ig_verts() const
    {
        return 5;
    }

    moris_index
    get_num_ig_cells() const
    {
        return 4;
    }

    moris_index
    get_num_verts_per_cell() const
    {
        return 3;
    }

    Matrix< IndexMat >
    get_new_node_on_parent_edge() const
    {
        return { {} };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_edge_edge_ordinal() const
    {
        return { {} };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_face() const
    {
        return {{}};
    }

    Matrix< IndexMat >
    get_new_node_on_parent_face_face_ordinal() const
    {
        return {{}};
    }

    Matrix< IndexMat >
    get_new_node_in_parent_cell() const
    {
        return { { 0 } };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_cell_cell_ordinal() const
    {
        return { { 0 } };
    }

    Vertex_Ancestry
    get_vertex_ancestry() const
    {
        return Vertex_Ancestry(
            { 0, 1, 2, 3, 0 },
            { EntityRank::NODE,
              EntityRank::NODE,
              EntityRank::NODE,
              EntityRank::NODE,
              EntityRank::ELEMENT } );
    }

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coordinates_wrt_parent() const
    {
        return { { { 0.0, 0.0 } } };
    }

    Cell< Cell< moris_index > >
    get_ig_cell_to_vertex_connectivity() const
    {
        return { 
            { { 0, 1, 4 } },
            { { 1, 2, 4 } },
            { { 2, 3, 4 } },
            { { 3, 0, 4 } } };
    }

    moris_index
    get_parametric_dimension() const
    {
        return 2;
    }

    moris_index
    get_signature() const
    {
        return 420;
    }

    enum CellTopology
    get_ig_cell_topology() const
    {
        return CellTopology::TRI3;
    }
};

// -------------------------------------------------------------------------

class Regular_Subdivision_24_TETS : public Regular_Subdivision_Template
{
  public:
    Regular_Subdivision_24_TETS() {}

    moris_index
    get_num_new_nodes() const
    {
        return 7;
    }

    moris_index
    get_total_ig_verts() const
    {
        return 15;
    }

    moris_index
    get_num_ig_cells() const
    {
        return 24;
    }

    moris_index
    get_num_verts_per_cell() const
    {
        return 4;
    }

    Matrix< IndexMat >
    get_new_node_on_parent_edge() const
    {
        return { {} };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_edge_edge_ordinal() const
    {
        return { {} };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_face() const
    {
        return { { 0, 1, 2, 3, 4, 5 } };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_face_face_ordinal() const
    {
        return { { 0, 1, 2, 3, 4, 5 } };
    }

    Matrix< IndexMat >
    get_new_node_in_parent_cell() const
    {
        return { { 6 } };
    }

    Matrix< IndexMat >
    get_new_node_on_parent_cell_cell_ordinal() const
    {
        return { { 0 } };
    }

    Vertex_Ancestry
    get_vertex_ancestry() const
    {
        return Vertex_Ancestry(
            { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 0 },
            { EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::NODE,
                EntityRank::FACE,
                EntityRank::FACE,
                EntityRank::FACE,
                EntityRank::FACE,
                EntityRank::FACE,
                EntityRank::FACE,
                EntityRank::ELEMENT } );
    }

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coordinates_wrt_parent() const
    {
        return { { { 0.0, -1.0, 0.0 } },
            { { 1.0, 0.0, 0.0 } },
            { { 0.0, 1.0, 0.0 } },
            { { -1.0, 0.0, 0.0 } },
            { { 0.0, 0.0, -1.0 } },
            { { 0.0, 0.0, 1.0 } },
            { { 0.0, 0.0, 0.0 } } };
    }

    Cell< Cell< moris_index > >
    get_ig_cell_to_vertex_connectivity() const
    {
        return { { { 0, 8, 1, 14 } },
            { { 1, 8, 5, 14 } },
            { { 4, 5, 8, 14 } },
            { { 0, 4, 8, 14 } },
            { { 1, 9, 2, 14 } },
            { { 2, 9, 6, 14 } },
            { { 5, 6, 9, 14 } },
            { { 1, 5, 9, 14 } },
            { { 2, 10, 3, 14 } },
            { { 2, 6, 10, 14 } },
            { { 6, 7, 10, 14 } },
            { { 3, 10, 7, 14 } },
            { { 0, 3, 11, 14 } },
            { { 3, 7, 11, 14 } },
            { { 4, 11, 7, 14 } },
            { { 0, 11, 4, 14 } },
            { { 0, 1, 12, 14 } },
            { { 1, 2, 12, 14 } },
            { { 2, 3, 12, 14 } },
            { { 0, 12, 3, 14 } },
            { { 4, 13, 5, 14 } },
            { { 5, 13, 6, 14 } },
            { { 6, 13, 7, 14 } },
            { { 4, 7, 13, 14 } } };
    }

    moris_index
    get_parametric_dimension() const
    {
        return 3;
    }

    moris_index
    get_signature() const
    {
        return 10000;
    }

    enum CellTopology
    get_ig_cell_topology() const
    {
        return CellTopology::TET4;
    }
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

class Generated_Regular_Subdivision_Template
{
  public:
    moris_index                 mNumNewNodes;
    moris_index                 mNumTotalIgVerts;
    moris_index                 mNumIgCells;
    Cell< Matrix< DDRMat > >    mParamCoords;
    Cell< Cell< moris_index > > mIgCellToVertOrd;
    Vertex_Ancestry             mNewVertexAncestry;
    Cell< moris_index >         mVertexHash;

  public:
    Generated_Regular_Subdivision_Template() {}

    ~Generated_Regular_Subdivision_Template() {}

    moris_index
    get_num_new_nodes() const
    {
        return mNumNewNodes;
    }

    moris_index
    get_total_ig_verts() const
    {
        return mNumTotalIgVerts;
    }

    moris_index
    get_num_ig_cells() const
    {
        return mNumIgCells;
    }

    const Vertex_Ancestry*
    get_new_vertex_ancestry() const
    {
        return &mNewVertexAncestry;
    }

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coordinates_wrt_parent() const
    {
        return mParamCoords;
    }

    Cell< Cell< moris_index > >
    get_ig_cell_to_vertex_connectivity() const
    {
        return mIgCellToVertOrd;
    }
};

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

class Regular_Subdivision_Interface : public Decomposition_Algorithm
{
  private:
    std::shared_ptr< Regular_Subdivision_Template >                          mRegularSubdivisionTemplate;
    moris::Cell< std::shared_ptr< Generated_Regular_Subdivision_Template > > mGeneratedTemplate;

    Integration_Mesh_Generation_Data* mMeshGenerationData;
    Decomposition_Data*               mDecompositionData;
    Cut_Integration_Mesh*             mCutIntegrationMesh;
    moris::mtk::Mesh*                 mBackgroundMesh;
    Integration_Mesh_Generator*       mGenerator;
    moris::uint                       mNumTotalCells = 0;


  public:

    Regular_Subdivision_Interface( ParameterList& aParameterList, enum CellTopology aCellTopology );

    ~Regular_Subdivision_Interface() {}

    bool has_geometric_independent_vertices() const;

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


    enum Decomposition_Algorithm_Type
    get_algorithm_type() const;

    // template functions

    moris_index
    get_num_new_nodes() const
    {
        return mRegularSubdivisionTemplate->get_num_new_nodes();
    }
    moris_index
    get_total_ig_verts() const
    {
        return mRegularSubdivisionTemplate->get_total_ig_verts();
    }
    moris_index
    get_num_ig_cells() const
    {
        return mRegularSubdivisionTemplate->get_num_ig_cells();
    }

    Matrix< IndexMat >
    get_new_node_on_parent_edge() const
    {
        return mRegularSubdivisionTemplate->get_new_node_on_parent_edge();
    }
    Matrix< IndexMat >
    get_new_node_on_parent_edge_edge_ordinal() const
    {
        return mRegularSubdivisionTemplate->get_new_node_on_parent_edge_edge_ordinal();
    }

    Matrix< IndexMat >
    get_new_node_on_parent_face() const
    {
        return mRegularSubdivisionTemplate->get_new_node_on_parent_face();
    }
    Matrix< IndexMat >
    get_new_node_on_parent_face_face_ordinal() const
    {
        return mRegularSubdivisionTemplate->get_new_node_on_parent_face_face_ordinal();
    }

    Matrix< IndexMat >
    get_new_node_in_parent_cell() const
    {
        return mRegularSubdivisionTemplate->get_new_node_in_parent_cell();
    }
    Matrix< IndexMat >
    get_new_node_on_parent_cell_cell_ordinal() const
    {
        return mRegularSubdivisionTemplate->get_new_node_on_parent_cell_cell_ordinal();
    };

    Cell< Matrix< DDRMat > >
    get_new_vertex_parametric_coordinates_wrt_parent() const
    {
        return mRegularSubdivisionTemplate->get_new_vertex_parametric_coordinates_wrt_parent();
    }
    Cell< Cell< moris_index > >
    get_ig_cell_to_vertex_connectivity() const
    {
        return mRegularSubdivisionTemplate->get_ig_cell_to_vertex_connectivity();
    }

    moris_index
    get_parametric_dimension() const
    {
        return mRegularSubdivisionTemplate->get_parametric_dimension();
    }

    moris_index
    get_signature() const
    {
        return mRegularSubdivisionTemplate->get_signature();
    }

    enum CellTopology
    get_ig_cell_topology() const
    {
        return mRegularSubdivisionTemplate->get_ig_cell_topology();
    }

  private:
    void
    make_new_vertex_requests(
        Child_Mesh_Experimental*            aChildMesh,
        Cut_Integration_Mesh*               aCutIntegrationMesh,
        Regular_Subdivision_Interface*      aRegularSubdivisionInterface,
        Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
        moris::mtk::Mesh*                   aBackgroundMesh,
        Decomposition_Data*                 aDecompositionData );
    void
    make_new_vertex_requests_trivial(
        Child_Mesh_Experimental*            aChildMesh,
        Cut_Integration_Mesh*               aCutIntegrationMesh,
        Regular_Subdivision_Interface*      aRegularSubdivisionInterface,
        Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
        moris::mtk::Mesh*                   aBackgroundMesh,
        Decomposition_Data*                 aDecompositionData );
    void
    make_new_vertex_requests_octree(
        Child_Mesh_Experimental*            aChildMesh,
        Cut_Integration_Mesh*               aCutIntegrationMesh,
        Regular_Subdivision_Interface*      aRegularSubdivisionInterface,
        Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
        moris::mtk::Mesh*                   aBackgroundMesh,
        Decomposition_Data*                 aDecompositionData );

    void
    generate_new_node_parent_information_ijk_mesh(
        Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
        moris::uint                         aNumIgCells,
        Child_Mesh_Experimental*            aChildMesh );
};


}// namespace xtk
#endif