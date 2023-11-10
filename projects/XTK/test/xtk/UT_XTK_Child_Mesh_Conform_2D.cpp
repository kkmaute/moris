/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Child_Mesh_Conform_2D.cpp
 *
 */

#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_plus.hpp"

#include "cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "fn_GEN_Triangle_Geometry.hpp" // For surface normals
#include "fn_equal_to.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

namespace xtk
{
class Permutations_2D
{
public:
    Permutations_2D()
    {
        load_all_2d_permutations();
    }

    size_t get_num_permutations()
    {
        return mNumPermutations;
    }

    Cell<size_t>
    get_permutation(size_t aIndex)
    {
        return mPermutations(aIndex);
    }

private:
    size_t mNumPermutations;
    Cell<Cell<size_t>> mPermutations;

    void
    load_all_2d_permutations()
    {
        mPermutations = {{0,1},{0,2},{1,2}};
        mNumPermutations = mPermutations.size();
    }
};

void
get_midside_coordinate(moris::moris_index const & aEdgeIndex,
                       moris::Matrix< moris::DDRMat > const & aNodeCoordinates,
                       moris::Matrix< moris::DDRMat >       & aMidEdgeCoordinate)
{
    if( aEdgeIndex == 0)
    {
        aMidEdgeCoordinate(0,0) = 0.5*(aNodeCoordinates(0,0) + aNodeCoordinates(1,0));
        aMidEdgeCoordinate(0,1) = 0.5*(aNodeCoordinates(0,1) + aNodeCoordinates(1,1));
    }

    else if( aEdgeIndex == 1)
    {
        aMidEdgeCoordinate(0,0) = 0.5*(aNodeCoordinates(1,0) + aNodeCoordinates(2,0));
        aMidEdgeCoordinate(0,1) = 0.5*(aNodeCoordinates(1,1) + aNodeCoordinates(2,1));
    }

    else if( aEdgeIndex == 2)
    {
        aMidEdgeCoordinate(0,0) = 0.5*(aNodeCoordinates(0,0) + aNodeCoordinates(2,0));
        aMidEdgeCoordinate(0,1) = 0.5*(aNodeCoordinates(0,1) + aNodeCoordinates(2,1));
    }
    else
    {
        std::cout<<"UNDEFINED EDGE"<<std::endl;
    }
}

void
setup_node_coordinates(moris::moris_index const & tEdge0,
                       moris::Matrix< moris::DDRMat > & aNodeCoordinates)
{
    aNodeCoordinates = moris::Matrix< moris::DDRMat >(4,2);

    aNodeCoordinates(0,0) = 0.0; aNodeCoordinates(0,1) = 0.0;
    aNodeCoordinates(1,0) = 1.0; aNodeCoordinates(1,1) = 0.0;
    aNodeCoordinates(2,0) = 1.0; aNodeCoordinates(2,1) = 1.0;

    moris::Matrix< moris::DDRMat > tEdgeNodeCoordinates(1,2,75);
    get_midside_coordinate(tEdge0,aNodeCoordinates, tEdgeNodeCoordinates);
    aNodeCoordinates(3,0) = tEdgeNodeCoordinates(0,0);     aNodeCoordinates(3,1) =  tEdgeNodeCoordinates(0,1);
}

void
setup_node_coordinates(moris::moris_index const & tEdge0,
                       moris::moris_index const & tEdge1,
                       moris::Matrix< moris::DDRMat > & aNodeCoordinates)
{
    aNodeCoordinates = moris::Matrix< moris::DDRMat >(5,2);

    aNodeCoordinates(0,0) = 0.0; aNodeCoordinates(0,1) = 0.0;
    aNodeCoordinates(1,0) = 1.0; aNodeCoordinates(1,1) = 0.0;
    aNodeCoordinates(2,0) = 1.0; aNodeCoordinates(2,1) = 1.0;

    moris::Matrix< moris::DDRMat > tEdgeNodeCoordinates(1,2,75);
    get_midside_coordinate(tEdge0,aNodeCoordinates, tEdgeNodeCoordinates);
    aNodeCoordinates(3,0) = tEdgeNodeCoordinates(0,0);     aNodeCoordinates(3,1) =  tEdgeNodeCoordinates(0,1);

    get_midside_coordinate(tEdge1, aNodeCoordinates, tEdgeNodeCoordinates);
    aNodeCoordinates(4,0) = tEdgeNodeCoordinates(0,0);     aNodeCoordinates(4,1) =  tEdgeNodeCoordinates(0,1);

}

TEST_CASE(" 2-D Conformal Decomposition Templates ","[2D_Conf_Temp]")
{
    // base cell information
    moris::Matrix< moris::IndexMat > tNodeIndex({{0,1,2}});
    moris::Matrix< moris::IdMat > tNodeIds({{1,2,3,4,5}});
    moris::Matrix< moris::IndexMat > tElementsAncestry({{0}}); // Not used
    moris::Matrix< moris::DDSTMat  > tElementNodeParentRanks(1,3,0);
    moris::Matrix< moris::IndexMat > tParentEdgeInds({{0,1,2}});
    moris::Matrix< moris::DDSTMat  > tParentEdgeRanks(1,3,1);

    Permutations_2D tPermutations;
    size_t tNumPermutations = tPermutations.get_num_permutations();
    for(size_t iPerm = 0; iPerm<tNumPermutations; iPerm++)
    {
        // Initialize Template
        Mesh_Modification_Template tMeshTemplate(tElementsAncestry(0,0),
                                                 0,
                                                 tNodeIndex,
                                                 tNodeIndex,
                                                 tElementNodeParentRanks,
                                                 tParentEdgeInds,
                                                 tParentEdgeRanks,
                                                 {{}},
                                                 {{}},
                                                 TemplateType::TRI_3);

        // Initialize child mesh with template (in this case a tet4)
        Child_Mesh tChildMesh(tMeshTemplate);

        tChildMesh.add_new_geometry_interface(0);

        // add new node indices
        tChildMesh.add_node_indices({{3,4}});
        tChildMesh.add_node_ids(tNodeIds);

        // select template
        Cell<size_t> tCurrentPermutation = tPermutations.get_permutation(iPerm);
        moris::moris_index tEdge0  = tCurrentPermutation(0);
        moris::moris_index tEdge1  = tCurrentPermutation(1);

        // Set up node coordinates
        moris::Matrix< moris::DDRMat > tNodeCoords;
        setup_node_coordinates(tEdge0,tEdge1,tNodeCoords);

        tChildMesh.add_entity_to_intersect_connectivity(3,tEdge0,true);
        tChildMesh.add_entity_to_intersect_connectivity(4,tEdge1,true);

        tChildMesh.modify_child_mesh(TemplateType::CONFORMAL_TRI3);

        moris_id tElemId = 1;
        tChildMesh.set_child_element_ids(tElemId);

        moris::Matrix< moris::IdMat > tElementIds  = tChildMesh.get_element_ids();
        moris::Matrix< moris::IdMat > tElementNode = tChildMesh.get_element_to_node_global();

    }
}

TEST_CASE(" 2-D Conformal Coincident Decomposition Templates ","[2D_Conf_Coin_Temp]")
{
    // base cell information
    moris::Matrix< moris::IndexMat > tNodeIndex({{0,1,2}});
    moris::Matrix< moris::IdMat > tNodeIds({{1,2,3,4}});
    moris::Matrix< moris::IndexMat > tElementsAncestry({{0}}); // Not used
    moris::Matrix< moris::DDSTMat  > tElementNodeParentRanks(1,3,0);
    moris::Matrix< moris::IndexMat > tParentEdgeInds({{0,1,2}});
    moris::Matrix< moris::DDSTMat  > tParentEdgeRanks(1,3,1);

    // Permutations_2D tPermutations;
    // size_t tNumPermutations = tPermutations.get_num_permutations();
    for(size_t iPerm = 0; iPerm<3; iPerm++)
    {
        // Initialize Template
        Mesh_Modification_Template tMeshTemplate(tElementsAncestry(0,0),
                                                 0,
                                                 tNodeIndex,
                                                 tNodeIndex,
                                                 tElementNodeParentRanks,
                                                 tParentEdgeInds,
                                                 tParentEdgeRanks,
                                                 {{}},
                                                 {{}},
                                                 TemplateType::TRI_3);

        // Initialize child mesh with template (in this case a tet4)
        Child_Mesh tChildMesh(tMeshTemplate);

        tChildMesh.add_new_geometry_interface(0);

        // add new node indices
        tChildMesh.add_node_indices({{3}});
        tChildMesh.add_node_ids(tNodeIds);

        // select template
        moris::moris_index tEdge0  = iPerm;

        // Set up node coordinates
        moris::Matrix< moris::DDRMat > tNodeCoords;
        setup_node_coordinates(tEdge0,tNodeCoords);

        tChildMesh.add_entity_to_intersect_connectivity(3,tEdge0,true);

        tChildMesh.modify_child_mesh(TemplateType::CONFORMAL_TRI3);

        moris_id tElemId = 1;
        tChildMesh.set_child_element_ids(tElemId);

        moris::Matrix< moris::IdMat > tElementIds  = tChildMesh.get_element_ids();
        moris::Matrix< moris::IdMat > tElementNode = tChildMesh.get_element_to_node_global();

    }
}
}

