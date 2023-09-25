/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Intersection_Detect.cpp
 *
 */

#include "cl_XTK_Intersection_Detect.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_create_edges_from_element_to_node.hpp"
#include "fn_generate_face_to_face.hpp"
#include "fn_unique.hpp"
#include<unordered_map>

namespace xtk
{

    Intersection_Detect::Intersection_Detect(
            std::shared_ptr<moris::mtk::Mesh_Manager> aMeshManager,
            Cut_Mesh & aCutMesh,
            moris::moris_index                   aMeshIndex,
            moris::ParameterList &        aParameterList):
                                                       mMeshManager(aMeshManager),
                                                       mCutMesh(aCutMesh),
                                                       mMeshIndex(aMeshIndex)
                                                       {
        //get the periodic mesh set names
        std::string tMeshSideSetNames = aParameterList.get<std::string>( "periodic_side_set_pair" );

        //store the names in mMeshSideSetPairs
        string_to_cell_of_cell(aParameterList.get<std::string>( "periodic_side_set_pair" ), mMeshSideSetPairs);
                                                       }

    //--------------------------
    void
    Intersection_Detect::find_intersections()
    {

        moris::mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh( mMeshIndex );

        moris::mtk::Set* tSet1 = tIntegrationMesh->get_set_by_name( "SideSet_2_c_p1" );

        // get clusters in the the first set
        moris::Cell< moris::mtk::Cluster const* > tSetClusters = tSet1->get_clusters_on_set();

        moris::Cell<moris::mtk::Cell const *> const & tCells = tSetClusters(0)->get_primary_cells_in_cluster();

        for( uint i = 0 ; i < tCells.size() ; i ++)
          {
              std::cout<<tCells(i)<<std::endl;
          }

        moris::Matrix<moris::IndexMat> tCellOrds = tSetClusters(0)->get_cell_side_ordinals();

        print(tCellOrds,"tCellOrdinals");

        moris::Matrix<moris::IdMat> ElementToNode( tCells.size(), 3 );

        for( uint i = 0 ; i < tCells.size() ; i++ )
        {
            moris::Cell<moris::mtk::Vertex const *> tVerticess = tCells(i)->get_vertices_on_side_ordinal( tCellOrds(i) );

            for ( uint j = 0 ; j < tVerticess.size() ; j++ )
            {
                print(tVerticess(j)->get_coords(), "Coordinates");
                ElementToNode( i, j ) = tVerticess(j)->get_id();
                moris::Matrix <moris::DDRMat > tCoord = tVerticess(j)->get_coords();
            }
        }

        print(ElementToNode, "ElementToNode");
//        create_edges_from_element_to_node(
//                    mtk::CellTopology                   aElementTopology,
//                    Integer                             aNumNodes,
//                    moris::Matrix< moris::IdMat > const & aElementToNode,
//                    moris::Matrix< moris::IdMat >       & aElementToEdge,
//                    moris::Matrix< moris::IdMat >       & aEdgeToNode,
//                    moris::Matrix< moris::IdMat >       & aNodeToEdge,
//                    moris::Matrix< moris::IdMat >       & aEdgeToElement)
        moris::Matrix< moris::IdMat >        aElementToEdge;
        moris::Matrix< moris::IdMat >        aEdgeToNode;
        moris::Matrix< moris::IdMat >        aNodeToEdge;
        moris::Matrix< moris::IdMat >        aEdgeToElement;

       // ElementToNode = {{0,1,2},{2,3,1}};

        moris::Matrix<moris::IdMat> tUniqueElementIds;
        unique(ElementToNode,tUniqueElementIds);

        print(tUniqueElementIds,"tUniqueElementIds");

        //construct a map
        std::unordered_map<moris_id,moris_index> tGlobaltoLocalmap;

        for(moris::uint iB =0; iB < tUniqueElementIds.numel() ; iB++)
        {

            tGlobaltoLocalmap[tUniqueElementIds(iB)] = iB;
        }

        moris::Matrix < moris::DDRMat > tCoords( 2, tUniqueElementIds.numel() );
        //assemble coordinates of the matrix
        for (uint i = 0 ; i < tUniqueElementIds.numel() ; i++ )
        {
            moris::mtk::Vertex const & tVertex = tIntegrationMesh->get_mtk_vertex( tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id(  tUniqueElementIds(i) , mtk::EntityRank::NODE ) );

            print(tVertex.get_coords(),"tVertex.get_coords()");

            tCoords.get_column(i)(0) = tVertex.get_coords()(1);
            tCoords.get_column(i)(1) = tVertex.get_coords()(2);
        }

        moris::Matrix<moris::IndexMat> tLocalElemToNode( ElementToNode.n_rows(), ElementToNode.n_cols() );
//        for( uint i = 0 ; i < ElementToNode.numel() ; i++ )
//        {
//            //auto tIter = tGlobaltoLocalmap.find( ElementToNode( i ) );
//
//            //std::cout<<tGlobaltoLocalmap.at( ElementToNode( i ) )<<std::endl;
//            std::cout<<"  "<<tGlobaltoLocalmap[tUniqueElementIds( i )]<<std::endl;
//
//            tLocalElemToNode ( i ) = tGlobaltoLocalmap.at( ElementToNode( i ) );
//        }

        for( uint tRowNum = 0 ; tRowNum < ElementToNode.n_rows() ; tRowNum++ )
        {
            for( uint tColNum = 0 ; tColNum < ElementToNode.n_cols() ; tColNum++ )
            {
            //auto tIter = tGlobaltoLocalmap.find( ElementToNode( i ) );

            std::cout<<tGlobaltoLocalmap.at( ElementToNode( tRowNum,tColNum ) )<<std::endl;
            //std::cout<<"  "<<tGlobaltoLocalmap[tUniqueElementIds( i )]<<std::endl;

            tLocalElemToNode ( tRowNum, tColNum ) = tGlobaltoLocalmap.at( ElementToNode( tRowNum, tColNum ) );
            }
        }

        print(tLocalElemToNode,"tLocalElemToNode");

        create_edges_from_element_to_node(
                    mtk::CellTopology::TRI3,
                    ElementToNode.n_cols()*ElementToNode.n_rows(),
                    tLocalElemToNode,
                    aElementToEdge,
                    aEdgeToNode,
                    aNodeToEdge,
                    aEdgeToElement);

        print(aNodeToEdge,"aNodeToEdge");
        print(aEdgeToElement, "aEdgeToElement");
        print(aElementToEdge, "aElementToEdge");

        moris::Matrix< moris::IndexMat > tElemToElemLocal = generate_face_to_face(aElementToEdge,ElementToNode.n_rows(),ElementToNode.n_cols(), std::numeric_limits<moris::moris_index>::max()) ;

        print(tElemToElemLocal ,"tElemToElemLocal") ;

        moris::moris_index aChildMeshIndex = 0.0;
        moris::moris_index aParentFaceIndex = 3.0 ;
        moris::Matrix<moris::IdMat> aChildrenElementId;
        moris::Matrix<moris::IndexMat> aChildrenElementCMInd;
        moris::Matrix<moris::IndexMat> aFaceOrdinal;

        std::cout<<"NumChildMesh: "<<mCutMesh.get_num_child_meshes()<<std::endl;

        for (uint i = 0 ; i < mCutMesh.get_num_child_meshes(); i ++ )
        {
            mCutMesh.get_child_elements_connected_to_parent_facet(
                    i,
                    aParentFaceIndex,
                    aChildrenElementId,
                    aChildrenElementCMInd,
                    aFaceOrdinal);

//            print(aChildrenElementId,"aChildrenElementId");
//            print(aChildrenElementCMInd, "aChildrenElementCMInd");
//            print(aFaceOrdinal, "aFaceOrdinal");
        }

        Child_Mesh const & CM0 = mCutMesh.get_child_mesh(aChildMeshIndex);

        moris::Cell<moris::mtk::Vertex const *> const & tvertices = CM0.get_vertices();

        std::cout<<tvertices.size()<<"tvertices.size()" << std::endl;

//        for( uint i = 0 ; i < tvertices.size() ; i++ )
//        {
//            print( tvertices(i)->get_coords(), "MatrixCoord");
//        }

    }
}

