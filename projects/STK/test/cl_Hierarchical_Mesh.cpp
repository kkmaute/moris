// Third-party header files.
#include <catch.hpp>
#include <iostream>
#include <mpi.h>
#include <string>

// MORIS project header files.
#include "algorithms.hpp"
#include "cl_Hierarchical_Mesh_Main.hpp" // STK/src/Hierarchical

#include "cl_Mesh.hpp" // MTK/src
#include "cl_Database.hpp" // MTK/src
#include "cl_STK_Implementation.hpp" // STK/src
#include "cl_HMR_Implementation.hpp" // STK/src/Hierarchical


using namespace moris;
// ----------------------------------------------------------------------------

TEST_CASE("HierarchicalMesh", "[moris],[mesh],[cl_Mesh],[Mesh]")
{

    SECTION( "Compute the T-Matrix for a specific element ")
                                                {
        if( par_rank() == 0)
        {
            if( par_size() == 1)
            {
                // This example creates a hierarchical mesh
                // It demonstrates the functionality to create a T-Matrix and the ID-Field
                uint tModelDim = 2; // ModelDimension
                uint tPolynomial = 1; // Polynomial degree of the shape functions
                Mat<uint> tNumElements ={ { 4, 4, 0 } }; // Number of Elements in each direction
                Mat<real> tDimensions ={ { 4.0, 4.0, 0.0 } }; // Number of Elements in each direction
                // Construct tHMR
                Hierarchical_Mesh_Main tHMR( tModelDim, tPolynomial, tNumElements, tDimensions );
                tHMR.create_mesh(); // Create a mesh (All data gets calculated to compute information or output a mesh)
                tHMR.activate_basisfunction(); // Activate basis functions for the current active elements (initial mesh on level zero is active)
                tHMR.give_Tmatrix_and_IdField( 28 ); // Create the T-Matrix and the ID-Field for element 28
                REQUIRE( tHMR.mElementData.IdField( 0 ) == 32 );
                REQUIRE( tHMR.mElementData.IdField( 1 ) == 33 );
                REQUIRE( tHMR.mElementData.IdField( 2 ) == 39 );
                REQUIRE( tHMR.mElementData.IdField( 3 ) == 40 );
                //See the T-Matrix and the ID-Field for Element 28
                //tHMR.mElementData.TMatrix.print("T-Matrix of basis functions");
                //tHMR.mElementData.IdField.print("ID");
                Mat<real> tSpecificPhysicalPoint = { { 1.5, 1.5, 0.0} };
                 // Nodal Interpolation for the physical coordinate: 1.5, 1.5, 0.0
                tHMR.give_Tmatrix_and_IdField_Specific_PhysicalCoord( tSpecificPhysicalPoint );
                REQUIRE( tHMR.mElementData.IdField( 0 ) == 16 );
                REQUIRE( tHMR.mElementData.IdField( 1 ) == 17 );
                REQUIRE( tHMR.mElementData.IdField( 2 ) == 23 );
                REQUIRE( tHMR.mElementData.IdField( 3 ) == 24 );
                REQUIRE( tHMR.mElementData.TMatrix( 0 ) == 0.25 );
                REQUIRE( tHMR.mElementData.TMatrix( 1 ) == 0.25 );
                REQUIRE( tHMR.mElementData.TMatrix( 2 ) == 0.25 );
                REQUIRE( tHMR.mElementData.TMatrix( 3 ) == 0.25 );
                Mat<real> tSpecificNaturalPoint = { { 0.5, 0.5, 0.0 } };
                 // Element Id 14 and the natural position: 0.5, 0.5, 0.0
                tHMR.give_Tmatrix_and_IdField_Specific_NaturalCoord( 14, tSpecificNaturalPoint );
                REQUIRE( tHMR.mElementData.IdField( 0 ) == 16 );
                REQUIRE( tHMR.mElementData.IdField( 1 ) == 17 );
                REQUIRE( tHMR.mElementData.IdField( 2 ) == 23 );
                REQUIRE( tHMR.mElementData.IdField( 3 ) == 24 );
                REQUIRE( tHMR.mElementData.TMatrix( 0 ) == 0.25 );
                REQUIRE( tHMR.mElementData.TMatrix( 1 ) == 0.25 );
                REQUIRE( tHMR.mElementData.TMatrix( 2 ) == 0.25 );
                REQUIRE( tHMR.mElementData.TMatrix( 3 ) == 0.25 );
            }
        }
                                                }

    SECTION( "Compute lagrange basis from position and vice versa ")
    {
        //This example works only in serial, with 2 or 4 processors (because of the decomposition of the mesh)
        // if( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        // At the moment Hierarchical Mesh does not work in parallel
        if( par_size() == 1 )
        {
            //Create a mesh from data
            uint tModelDim = 2;
            uint tPolynomial = 2;
            Mat<uint> tNumberOfElementsPerDirection = { { 2 }, { 2 }, { 0 } };
            Mat<real> tDomainDimensions = { { 3.0, 3.0, 3.0 } };
            //Create object with data
            Hierarchical_Mesh_Main tHMR( tModelDim, tPolynomial, tNumberOfElementsPerDirection, tDomainDimensions );
            //Compute tBasisId from level zero and position i = 0, j = 0
            uint tBasis = tHMR.give_lagrange_basis_of_position( 0, { { 0 }, { 0 } } );
            REQUIRE( tBasis == 0 );
            //Compute tBasisId from level zero and position i = 0, j = 0
            tBasis = tHMR.give_lagrange_basis_of_position( 0, { { 4 }, { 4 } } );
            REQUIRE( tBasis == 24 );
            //Compute the position with the basis ID
            Mat<uint> tPositionOfBasis = tHMR.give_position_of_lagrange_basis( 0 );
            REQUIRE( tPositionOfBasis( 0 ) == 0 );
            REQUIRE( tPositionOfBasis( 1 ) == 0 );
            //Compute the position with the basis ID
            tPositionOfBasis = tHMR.give_position_of_lagrange_basis( 24 );
            REQUIRE( tPositionOfBasis( 0 ) == 4 );
            REQUIRE( tPositionOfBasis( 1 ) == 4 );
        }
    }

    SECTION( "Create a mesh and dump it into MTK")
    {
        //This example works only in serial, with 2 or 4 processors (because of the decomposition of the mesh)
        //if( par_size() == 1 || par_size() == 2 || par_size() == 4 )

        // At the moment Hierarchical Mesh does not work in parallel
        if( par_size() == 1 )
        {
            uint tModelDim = 2; // ModelDimension
            uint tPolynomial = 1; // Polynomial degree of the shape functions
            Mat<uint> tNumElements ={ { 2, 2, 0 } }; // Number of Elements in each direction
            Mat<real> tDimensions ={ { 4.0, 4.0, 0.0 } }; // Number of Elements in each direction
            Hierarchical_Mesh_Main tHMR(tModelDim,tPolynomial,tNumElements,tDimensions); // Construct tHMR
            //Create a side set
            Mat<uint> tSideSet = { { 1 }, { 2 }, { 3 }, { 4 }, { 0 }, { 0 } };
            //Set side set for HMR
            tHMR.set_side_set( tSideSet );
            tHMR.create_mesh(); // Create a mesh (All data gets calculated to compute information or output a mesh)
            //Deactivate the arbitrary element in the middle of these 3 elements
            Mat<uint> tDeac = { { 6 } };
            tHMR.mElementData.DeactivateElement = tDeac;
            //Refine, based on the vector of elements, which needs to be refined
            tHMR.hierarchical_element_refinement();
            // Activate basis functions for the current active elements (initial mesh on level zero is active)
            tHMR.activate_basisfunction();
            //Create mesh data for MTK
            tHMR.create_mesh_data_new();
            //Create MTK file
            tHMR.create_MTK_file("test_mesh.exo");
        }
    }

    SECTION( "Create a mesh, link it to MTK and ask questions to the mesh class")
    {
        // //This example works only in serial, with 2 or 4 processors (because of the decomposition of the mesh)
       //  if( par_size() == 1 || par_size() == 2 || par_size() == 4 )

        // At the moment Hierarchical Mesh does not work in parallel
        if( par_size() == 1 )
        {
            uint tModelDim = 2; // ModelDimension
            uint tPolynomial = 1; // Polynomial degree of the shape functions
            Mat<uint> tNumElements ={ { 2, 2, 0 } }; // Number of Elements in each direction
            Mat<real> tDimensions ={ { 4.0, 4.0, 0.0 } }; // Number of Elements in each direction
            // Construct tHMR
            Hierarchical_Mesh_Main tHMR( tModelDim, tPolynomial, tNumElements, tDimensions );
            //Create a side set
            Mat<uint> tSideSet = { { 1 }, { 2 }, { 3 }, { 4 }, { 0 }, { 0 } };
            //Set side set for HMR
            tHMR.set_side_set( tSideSet );
            tHMR.create_mesh(); // Create a mesh (All data gets calculated to compute information or output a mesh)
            // Activate basis functions for the current active elements (initial mesh on level zero is active)
            tHMR.activate_basisfunction();
            //Create mesh data for MTK
            tHMR.create_mesh_data_new();
            //Create tMesh from HMR
            HMR_Implementation* tHMRimplementation  = new HMR_Implementation( tHMR );
            mesh tMesh( tHMRimplementation );
            //Now you can use all the commands from cl_mesh
            //Get the number of elements
            uint tNumberOfElements = tMesh.get_num_elems();
            REQUIRE( tNumberOfElements == 4 );
        }
    }




    //    SECTION( "Functionality of Elements and Basis Functions in a Hierarchical Mesh ")
    //    {
    //        // This example provides all the functionality for a hierarhical mesh
    //        // It demonstrates the possibility to determine elements or basic functions from both sides in a tensorial grid.
    //        uint tModelDim = 2; // ModelDimension
    //        uint tPolynomial = 2; // Polynomial degree of the shape functions
    //        Mat<uint> tNumElements ={{4,4,4}}; // Number of Elements in each direction
    //        Mat<real> tDimensions ={{4.0,4.0,4.0}}; // Number of Elements in each direction
    //        Hierarchical_Mesh_Main mHMR(tModelDim,tPolynomial,tNumElements,tDimensions);
    //
    //        // Functionalities for Elements
    //        Mat<uint> tElementNeighbor = mHMR.give_neighbor_of_element(aElement,aBuffer); // Give neighbors of an element
    //        uint tElementPosition = mHMR.give_element_of_position(aLevel,aPosition);
    //        Mat<uint> tPositionOfElement = mHMR.give_position_of_element(aElement);
    //        uint tElementLevel = mHMR.give_element_level(aElement);
    //        uint tParent = mHMR.give_parent_of_element(aElement);
    //        uint tParentX = mHMR.give_parent_of_level_x(aElement,aLevel);
    //        uint tNumElem = mHMR.give_number_of_elements(aLevel);
    //        Mat<uint> tChild = mHMR.give_children_of_element(aElement);
    //        Mat<uint> tBasisOfElement = mHMR.give_basis_of_element(aElement);
    //        uint tElementOwner = mHMR.give_element_owner(aElement);
    //        Mat<uint> Elementshare = mHMR.give_element_share(aElement); // Creates a list of sharing elements with Aura
    //        //
    //        //        // Functionalities for Basis
    //        //        uint tNumBasis = mHMR.give_number_of_basis(aLevel);
    //        //        Mat<uint> tNeighBasis = mHMR.give_neighbor_of_basis(aBasis,aBuffer);
    //        //        uint tBasLevel = mHMR.give_basis_level(aBasis);
    //        //        Mat<uint> tPositionOfBasis = mHMR.give_position_of_basis(aBasis);
    //        //        uint tBasisOfPosition = mHMR.give_basis_of_position(aLevel,aPosition);
    //        //        Mat<uint> tElementOfBasis = mHMR.give_element_of_basis(aBasis);
    //        //
    //        //        // Functionalities for Edges
    //        //        uint tNumEdgex = mHMR.give_number_of_edges_x(aLevel);
    //        //        uint tNumEdgey = mHMR.give_number_of_edges_y(aLevel);
    //        //        uint tNumEdgez = mHMR.give_number_of_edges_z(aLevel);
    //        //        uint tNumEdge = mHMR.give_number_of_edges(aLevel);
    //        //        uint tEdgexofPos = mHMR.give_edge_x_of_position(aLevel,aPosition);
    //        //        uint tEdgeyofPos = mHMR.give_edge_y_of_position(aLevel,aPosition);
    //        //        uint tEdgezofPos = mHMR.give_edge_z_of_position(aLevel,aPosition);
    //        //        uint tEdgeLEvel = mHMR.give_edge_level(aEdge);
    //        //        Mat<uint> tEdgesofElement = mHMR.give_element_edges(aElement);
    //        //        Mat<uint> tEdgePosition = mHMR.give_edge_position(aEdge);
    //        //        Mat<uint> tEdgeNodes = mHMR.give_edge_nodes(aEdge);
    //        //        uint tEdgeOwner = mHMR.give_edge_owner(aEdge);
    //        //        Mat<uint> Edgeshare = mHMR.give_edge_share(aEdge);
    //        //
    //        //        // Functionalities for Faces
    //        //        uint tNumFacex = mHMR.give_number_of_faces_x(aLevel);
    //        //        uint tNumFacey = mHMR.give_number_of_faces_y(aLevel);
    //        //        uint tNumFacez = mHMR.give_number_of_faces_z(aLevel);
    //        //        uint tNumFace = mHMR.give_number_of_faces(aLevel);
    //        //        uint tFaceLevel = mHMR.give_face_level(aFace);
    //        //        uint tFacexofPos = mHMR.give_face_x_of_position(aLevel,aPosition);
    //        //        uint tFaceyofPos = mHMR.give_face_y_of_position(aLevel,aPosition);
    //        //        uint tFacezofPos = mHMR.give_face_z_of_position(aLevel,aPosition);
    //        //        Mat<uint> tFaceofElement = mHMR.give_element_faces(aElement);
    //        //        Mat<uint> tFacePos = mHMR.give_face_position(aFace);
    //        //        Mat<uint> tNodesofFac = mHMR.give_face_nodes(aFace);
    //        //        uint tFaceOwner = mHMR.give_face_owner(aFace);
    //        //        Mat<uint> Faceshare = mHMR.give_face_share(aFace);
    //
    //        //Create a mesh in serial or parallel (Creates not a real mesh, but creates information to ask questions to the mesh)
    //        //        mHMR.create_mesh();
    //        //
    //        //        Mat<uint> tElementList = mHMR.give_element_on_proc();
    //        //        tBasisOfElement = mHMR.give_basis_of_element(tElementList(0));
    //        //        Mat<real> tCoordOfBasis = mHMR.give_coordinate_from_basis(tBasisOfElement(3));
    //
    //    }

}

//namespace moris
//{
//class HMR_XTK_Interface
//{
//public:
//    HMR_XTK_Interface(uint const & aModelDim,
//                      uint const & aPolynomial,
//                      Mat<uint> const & aNumElements,
//                      Mat<real> const & aModelDomain):
//                          mHMRMesh(aModelDim,aPolynomial,aNumElements,aModelDomain)
//    {
////        mHMRMesh.create_mesh();
//    }
//
//    Mat<uint>
//    get_nodes_connected_to_element(uint & aElementId)
//    {
//        Mat<uint> tElementBasis = mHMRMesh.give_basis_of_element(aElementId);
//        tensor_grid_to_topology_ordinals_nodes(tElementBasis);
//        return tElementBasis;
//    }
//
//    Mat<uint>
//    get_edges_connected_to_element(uint & aElementId)
//    {
//        Mat<uint> tElementEdgesTG = mHMRMesh.give_element_edges(aElementId);
//        Mat<uint> tElementEdgesTO = tensor_grid_to_topology_ordinals_edges(tElementEdgesTG);
//        return tElementEdgesTO;
//    }
//
//    Mat<uint>
//    get_faces_connected_to_element(uint & aElementId)
//    {
//        Mat<uint> tElementFacesTG = mHMRMesh.give_element_faces(aElementId);
//        Mat<uint> tElementFacesTO = tensor_grid_to_topology_ordinals_faces(tElementFacesTG);
//        return tElementFacesTO;
//    }
//
//private:
//    moris::Hierarchical_Mesh_Main mHMRMesh;
//
//
//    /**
//     * @param[in] aEntityRank - The rank of the entity to convert to a topology ordinal format
//     * @param[in] aConnectivity - Tensor format basis of a given entity
//     * @param[out] aConnectivity - Topology Ordinal Format of the given entity
//     */
//    void
//    tensor_grid_to_topology_ordinals_nodes(Mat<uint> & aConnectivity)
//    {
//        if(aConnectivity.length() == 8)
//        {
//            uint tVal = aConnectivity(2);
//            aConnectivity(2) = aConnectivity(3);
//            aConnectivity(3) = tVal;
//
//            tVal = aConnectivity(6);
//            aConnectivity(6) = aConnectivity(7);
//            aConnectivity(7) = tVal;
//        }
//
//        else if(aConnectivity.length() == 4)
//        {
//            uint tVal = aConnectivity(2);
//            aConnectivity(2) = aConnectivity(3);
//            aConnectivity(3) = tVal;
//        }
//        else
//        {
//            MORIS_LOG_ERROR<<"Invalid Connectivity Provided";
//        }
//    }
//
//    Mat<uint>
//    tensor_grid_to_topology_ordinals_edges(Mat<uint> & aInputConn)
//    {
//        Mat<uint> tOutputConn(aInputConn.length(),1);
//
//        if(aInputConn.length() == 12)
//        {
//            tOutputConn(0)  = aInputConn(0);
//            tOutputConn(1)  = aInputConn(5);
//            tOutputConn(2)  = aInputConn(1);
//            tOutputConn(3)  = aInputConn(4);
//            tOutputConn(4)  = aInputConn(2);
//            tOutputConn(5)  = aInputConn(7);
//            tOutputConn(6)  = aInputConn(3);
//            tOutputConn(7)  = aInputConn(6);
//            tOutputConn(8)  = aInputConn(8);
//            tOutputConn(9)  = aInputConn(9);
//            tOutputConn(10) = aInputConn(11);
//            tOutputConn(11) = aInputConn(10);
//        }
//
//        else if(aInputConn.length() == 4)
//        {
//            MORIS_LOG_ERROR<<"Not Implemented for a 2D Element yet";
//        }
//        else
//        {
//            MORIS_LOG_ERROR<<"Invalid Connectivity Provided";
//        }
//        return tOutputConn;
//    }
//
//    Mat<uint>
//    tensor_grid_to_topology_ordinals_faces(Mat<uint> & aInputConn)
//     {
//         Mat<uint> tOutputConn(aInputConn.length(),1);
//
//         if(aInputConn.length() == 6)
//         {
//             tOutputConn(0)  = aInputConn(3);
//             tOutputConn(1)  = aInputConn(2);
//             tOutputConn(2)  = aInputConn(4);
//             tOutputConn(3)  = aInputConn(1);
//             tOutputConn(4)  = aInputConn(5);
//             tOutputConn(5)  = aInputConn(6);
//         }
//         else
//         {
//             MORIS_LOG_ERROR<<"Invalid Connectivity Provided";
//         }
//         return tOutputConn;
//     }
//
//};
//}

//#include "cl_XTK_Model.hpp" // XTK/src/

//TEST_CASE("XTK to HMR Interface","[XTK_HMR]")
//{
//    moris::uint tDim=3; // Select the two or three dimensional case
//    moris::uint tPolynomial=1; // Select an arbitrary polynomial degree of the spline
//    moris::Mat<moris::uint> tNumElements = {{1},{1},{1}}; // Select the number of elements in each direction, x,y,z
//    moris::Mat<moris::real> tModelDomain = {{1},{1},{1}};
//
//
//    // Directly Creating the HMR
//    moris::Hierarchical_Mesh_Main tHierarchicalMesh(tDim,tPolynomial,tNumElements,tModelDomain);
////    tHierarchicalMesh.create_mesh();
//
//    // Wrapper around HMR
//    HMR_XTK_Interface tHMRWrap(tDim,tPolynomial,tNumElements,tModelDomain);
//
//
//
//    // STK Comparison
//    const std::string tFileName = "generated:1x1x1";
//    // Create MORIS mesh using MTK database
//    moris::mesh tMeshSTK( MeshType::MTK, tFileName );
//    moris::Mat<moris::uint> tElementNodeSTK = tMeshSTK.get_nodes_connected_to_element(1) ;
//    moris::Mat<moris::uint> tElementEdgeSTK = tMeshSTK.get_edges_connected_to_element(1) ;
//    moris::Mat<moris::uint> tElementFaceSTK = tMeshSTK.get_faces_connected_to_element(1) ;
//
//    for(moris::uint iEdge = 0; iEdge<tElementEdgeSTK.length(); iEdge++)
//    {
//        moris::Mat<moris::uint> tEdgeNode = tMeshSTK.get_nodes_connected_to_edge(tElementEdgeSTK(iEdge));
//        std::cout<<"Edge:" << iEdge<<" | ";
//        for(moris::uint iP = 0; iP<2; iP++)
//        {
//            std::cout<< tEdgeNode(iP) << " ";
//        }
//        std::cout<<std::endl;
//    }
//
//    for(moris::uint iFace = 0; iFace<tElementFaceSTK.length(); iFace++)
//    {
//        moris::Mat<moris::uint> tFaceNode = tMeshSTK.get_nodes_connected_to_face(tElementFaceSTK(iFace));
//        std::cout<<"Face:" << iFace<<" | ";
//        for(moris::uint iP = 0; iP<4; iP++)
//        {
//            std::cout<< tFaceNode(iP) << " ";
//        }
//        std::cout<<std::endl;
//    }
//
//
//    moris::uint tElementId = 0;
//    moris::uint tElementLevel = 0;
//
//
//    Mat<uint>tElementNodeWrap = tHMRWrap.get_nodes_connected_to_element(tElementId);
//    Mat<uint>tElementEdgeWrap = tHMRWrap.get_edges_connected_to_element(tElementId);
//    Mat<uint>tElementFaceWrap = tHMRWrap.get_faces_connected_to_element(tElementId);
//
//
//    tElementNodeWrap.print("Reordered Nodes");
//    moris::Mat<moris::uint> tElementNodeHMR = tHierarchicalMesh.give_basis_of_element(tElementId);
//    tElementNodeHMR.print("Element 1 Nodes HMR");
//    tElementNodeSTK.print("Element 1 Nodes STK");
//
//    moris::Mat<moris::uint> tElementEdgeHMR = tHierarchicalMesh.give_element_edges(tElementId);
//    tElementEdgeWrap.print("Reorded Edges");
//    tElementEdgeHMR.print("Element 1 Edges HMR");
//    tElementEdgeSTK.print("Element 1 Edges STK");
//
//    moris::Mat<moris::uint> tElementFaceHMR = tHierarchicalMesh.give_element_faces(tElementId);
//    tElementFaceWrap.print("Reordered Faces");
//    tElementFaceHMR.print("Element 1 Face HMR");
//    tElementFaceSTK.print("Element 1 Face STK");
//
//    moris::Mat<moris::uint> tFaceNodes = tHierarchicalMesh.give_face_nodes(tElementFaceHMR(0));
//    tFaceNodes.print("Face Nodes");
//
//    std::cout<<"Exit code 0"<<std::endl;
//
//

//}
