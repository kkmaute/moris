#ifndef MORIS_CL_XTK_REGULAR_SUBDIVISION_INTERFACE_HPP_
#define MORIS_CL_XTK_REGULAR_SUBDIVISION_INTERFACE_HPP_

#include "cl_XTK_Decomposition_Algorithm.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh;
    }
}

namespace xtk
{
struct Regular_Subdivision_Interface_Data
{
    public:
    Matrix<IndexMat> mNewNodesOnEdges;    
    Matrix<IndexMat> mNewNodesOnFaces;    
    Matrix<IndexMat> mNewNodesOnCells;    
    Matrix<IndexMat> mNewNodesOnEdgesOrd; 
    Matrix<IndexMat> mNewNodesOnFacesOrd; 
    Matrix<IndexMat> mNewNodesOnCellsOrd; 
    
    Cell<Cell<moris_index>> mNewCellToNodeConnectivity;

    Cell<Matrix<DDRMat>> mNewNodeXi;
    Matrix<DDRMat> mNXi;
    
};

class Integration_Mesh_Generator;
class Integration_Mesh_Generation_Data;
class Decomposition_Data;
class Cut_Integration_Mesh;


class Regular_Subdivision_24_TETS
{
    public:
    Regular_Subdivision_24_TETS(){}
    ~Regular_Subdivision_24_TETS(){}

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

    Matrix<IndexMat> get_new_node_on_parent_edge() const 
    {
        return {{}};
    }

    Matrix<IndexMat> 
    get_new_node_on_parent_edge_edge_ordinal() const
    {
        return {{}};
    }

    Matrix<IndexMat> get_new_node_on_parent_face() const
    {
        return {{0,1,2,3,4,5}};
    }

    Matrix<IndexMat> get_new_node_on_parent_face_face_ordinal() const 
    {
        return {{0,1,2,3,4,5}};
    }

    Matrix<IndexMat> get_new_node_in_parent_cell() const
    {
        return {{6}};
    }

    Matrix<IndexMat> 
    get_new_node_on_parent_cell_cell_ordinal() const 
    {
        return {{0}};
    }

    Cell<Matrix<DDRMat>> get_new_vertex_parametric_coordinates_wrt_parent() const
    {
        return {{{ 0.0, -1.0,  0.0}},
               {{ 1.0,  0.0,  0.0}},
               {{ 0.0,  1.0,  0.0}},
               {{-1.0,  0.0,  0.0}},
               {{ 0.0,  0.0, -1.0}},
               {{ 0.0,  0.0,  1.0}},
               {{ 0.0,  0.0,  0.0}}};
    }

    Cell<Cell<moris_index>>
    get_ig_cell_to_vertex_connectivity() const
    {
        return {{{0, 8, 1,  14}},
                {{1, 8, 5,  14}},
                {{4, 5, 8,  14}},
                {{0, 4, 8,  14}},
                {{1, 9, 2,  14}},
                {{2, 9, 6,  14}},
                {{5, 6, 9,  14}},
                {{1, 5, 9,  14}},
                {{2, 10, 3, 14}},
                {{2, 6, 10, 14}},
                {{6, 7, 10, 14}},
                {{3, 10, 7, 14}},
                {{0, 3, 11, 14}},
                {{3, 7, 11, 14}},
                {{4, 11, 7, 14}},
                {{0, 11, 4, 14}},
                {{0, 1, 12, 14}},
                {{1, 2, 12, 14}},
                {{2, 3, 12, 14}},
                {{0, 12, 3, 14}},
                {{4, 13, 5, 14}},
                {{5, 13, 6, 14}},
                {{6, 13, 7, 14}},
                {{4, 7, 13, 14}}};
    }

    moris_index 
    get_parametric_dimension() const
    {
        return 3;
    }

    moris_index get_signature() const 
    {
        return 10000;
    }

    enum CellTopology get_ig_cell_topology() const
    {
        return CellTopology::TET4;
    }
};


class Regular_Subdivision_Interface : public Decomposition_Algorithm
{
    private: 
    Regular_Subdivision_24_TETS mRegularSubdivisionTemplate;
    public:
    Regular_Subdivision_Interface(){}
    ~Regular_Subdivision_Interface(){}

    void 
    perform_impl_vertex_requests(Integration_Mesh_Generation_Data*  aMeshGenerationData,
                                 Decomposition_Data*                 aDecompositionData,
                                 Cut_Integration_Mesh*               aCutIntegrationMesh,
                                 moris::mtk::Mesh*                   aBackgroundMesh,
                                 Integration_Mesh_Generator*         aMeshGenerator);

    void
    perform_impl_generate_mesh(Integration_Mesh_Generation_Data* aMeshGenerationData,
                               Decomposition_Data*               aDecompositionData,
                               Cut_Integration_Mesh*             aCutIntegrationMesh,
                               moris::mtk::Mesh*                 aBackgroundMesh,
                               Integration_Mesh_Generator*       aMeshGenerator);

    

    enum Decomposition_Algorithm_Type
    get_algorithm_type() const;

    // template functions

    moris_index get_num_new_nodes() const { return mRegularSubdivisionTemplate.get_num_new_nodes() ;}
    moris_index get_total_ig_verts() const { return mRegularSubdivisionTemplate.get_total_ig_verts() ;}
    moris_index get_num_ig_cells() const { return mRegularSubdivisionTemplate.get_num_ig_cells() ;}

    Matrix<IndexMat> get_new_node_on_parent_edge() const { return mRegularSubdivisionTemplate.get_new_node_on_parent_edge() ;}
    Matrix<IndexMat> get_new_node_on_parent_edge_edge_ordinal() const  { return mRegularSubdivisionTemplate.get_new_node_on_parent_edge_edge_ordinal() ;}

    Matrix<IndexMat> get_new_node_on_parent_face() const  { return mRegularSubdivisionTemplate.get_new_node_on_parent_face() ;}
    Matrix<IndexMat> get_new_node_on_parent_face_face_ordinal() const  { return mRegularSubdivisionTemplate.get_new_node_on_parent_face_face_ordinal() ;}

    Matrix<IndexMat> get_new_node_in_parent_cell() const  { return mRegularSubdivisionTemplate.get_new_node_in_parent_cell() ;}
    Matrix<IndexMat> get_new_node_on_parent_cell_cell_ordinal() const  { return mRegularSubdivisionTemplate.get_new_node_on_parent_cell_cell_ordinal() ;};

    Cell<Matrix<DDRMat>> get_new_vertex_parametric_coordinates_wrt_parent() const { return mRegularSubdivisionTemplate.get_new_vertex_parametric_coordinates_wrt_parent() ;}
    Cell<Cell<moris_index>>  get_ig_cell_to_vertex_connectivity() const { return mRegularSubdivisionTemplate.get_ig_cell_to_vertex_connectivity() ;}

    moris_index get_parametric_dimension() const  { return mRegularSubdivisionTemplate.get_parametric_dimension() ;}

    moris_index get_signature() const { return mRegularSubdivisionTemplate.get_signature() ;}

    enum CellTopology get_ig_cell_topology() const { return mRegularSubdivisionTemplate.get_ig_cell_topology() ;}

    private:
    void
    make_new_vertex_requests(Child_Mesh_Experimental* aChildMesh,
                             Regular_Subdivision_Interface* aRegularSubdivisionInterface,
                             Regular_Subdivision_Interface_Data* aRegularSubdivisionInterfaceData,
                             moris::mtk::Mesh*   aBackgroundMesh,
                             Decomposition_Data* aDecompositionData);
};


}
#endif