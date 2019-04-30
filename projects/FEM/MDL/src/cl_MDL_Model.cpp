
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"                    //MTK/src

#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"               //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "../../INT/src/cl_FEM_Element_Bulk.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_Element_Block.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_print.hpp" // for check

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        Model::Model(       mtk::Mesh_Manager*            aMeshManager,
                      const uint                          aBSplineOrder,
                            Cell< Cell< fem::IWG_Type > > aIWGTypeList,
                            Cell< moris_index >           aSidesetList,
                            Cell< fem::BC_Type >          aSidesetBCTypeList  ) : mMeshManager( aMeshManager )
        {
            // start timer
            tic tTimer1;

            mDofOrder = aBSplineOrder;
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: initialize
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Get pointers to interpolation and integration mesh
            moris::moris_index tMeshPairIndex = 0;
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh = nullptr;
            mMeshManager->get_mesh_pair( tMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // ask mesh about number of nodes on proc
            luint tNumberOfNodes = tInterpolationMesh->get_num_nodes();

            // create node objects
            mNodes.resize(  tNumberOfNodes, nullptr );

            for( luint k = 0; k<tNumberOfNodes; ++k )
            {
                mNodes( k ) = new fem::Node( &tInterpolationMesh->get_mtk_vertex( k ) );
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created %u FEM nodes in %5.3f seconds.\n\n",
                        ( unsigned int ) tNumberOfNodes,
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1.5: create IWGs
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // number of IWGs to be created
            uint tNumOfIWGs = aIWGTypeList.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs for the problem considered
            mIWGs.resize( tNumOfIWGs );

            // loop over the IWG types
            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                mIWGs( i ).resize( aIWGTypeList( i ).size(), nullptr );

                for( uint Ki = 0; Ki < aIWGTypeList( i ).size(); Ki++)
                {
                    // create an IWG with the factory for the ith IWG type
                    mIWGs( i )( Ki ) = tIWGFactory.create_IWGs( aIWGTypeList( i )( Ki ) );
                }
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer2;

            // a factory to create the elements
            fem::Element_Factory tElementFactory;

            // get the number of element to create
            luint tNumberOfEquationObjects = tIntegrationMesh->get_num_elems()
                                           + tIntegrationMesh->get_sidesets_num_faces( aSidesetList );

            uint tNumberElementBlocks = 1 + aSidesetList.size();

            // create equation objects
            mElements.reserve( tNumberOfEquationObjects );

            mElementBlocks.resize( tNumberElementBlocks, nullptr );

            //  Create Blockset Elements ----------------------------------------
            std::cout<<" Create Blockset Elements "<<std::endl;
            //------------------------------------------------------------------------------

            // ask mesh about number of elements on proc
            moris::Cell<std::string> tBlockSetsNames = tIntegrationMesh->get_set_names( EntityRank::ELEMENT );

            // get cells in blockset (this needs to stay in scope somehow) // FIXME BLOCK OF CLUSTERS
            moris::Cell<mtk::Cell const*> tBlockSetCells = tIntegrationMesh->get_block_set_cells( tBlockSetsNames( 0 ) );

            // create new fem element block
            moris::uint tElementBlockCounter = 0;
            mElementBlocks( tElementBlockCounter ) = new fem::Element_Block( tBlockSetCells, fem::Element_Type::BULK, mIWGs( 0 ), mNodes );

            mElements.append( mElementBlocks( tElementBlockCounter++ )->get_equation_object_list() );

            //  Create Blockset Elements ----------------------------------------
            std::cout<<" Create Sideset Elements "<<std::endl;
            //------------------------------------------------------------------------------

            moris::Cell<std::string> tSideSetsNames = tIntegrationMesh->get_set_names( EntityRank::FACE );

            print(tSideSetsNames,"tSideSetsNames");

            for( luint Ik = 0; Ik < aSidesetList.size(); ++Ik )
            {
                // get the treated sideset name
                std::string tTreatedSidesetName = tSideSetsNames( aSidesetList( Ik ) );

                // create a cell of sideset element
                moris::Cell< mtk::Cell const * > tSideSetElement;

                // create a list of sideset ordinal
                Matrix< IndexMat >  aSidesetOrdinals;

                // get the treated sideset elements and ordinals
                tInterpolationMesh->get_sideset_cells_and_ords( tTreatedSidesetName, tSideSetElement, aSidesetOrdinals );

                mElementBlocks( tElementBlockCounter ) = new fem::Element_Block( tSideSetElement, fem::Element_Type::SIDESET, mIWGs( Ik + 1 ), mNodes );

                mElements.append( mElementBlocks( tElementBlockCounter++ )->get_equation_object_list() );

            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created %u FEM elements in %5.3f seconds.\n\n",
                        ( unsigned int ) tNumberOfEquationObjects,
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create Model Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer3;

            //--------------------------FIXME------------------------------------
            // This part should not be needed anymore when MTK has all the functionalities
            Matrix< IdMat > tCommTable;
            moris::map< moris::moris_id, moris::moris_index > tIdToIndMap;
            moris::uint tMaxNumAdofs;
            moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1 );

            if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
            {
                if ( mDofOrder == 0 )
                {
                    mDofOrder  = this->get_lagrange_order_from_mesh();
                }

                // get map from mesh
                tInterpolationMesh->get_adof_map( mDofOrder, mCoefficientsMap );

                tCommTable   = tInterpolationMesh->get_communication_table();
                tIdToIndMap  = mCoefficientsMap;
                tMaxNumAdofs = tInterpolationMesh->get_num_coeffs( mDofOrder );
                tDofTypes1( 0 ) = MSI::Dof_Type::L2;
            }
            else
            {
                tCommTable.set_size( 1, 1, 0 );
                tMaxNumAdofs = 1000000;
                tDofTypes1( 0 ) = MSI::Dof_Type::TEMP;
            }
            //--------------------------END FIXME--------------------------------

            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( mElementBlocks,
                                                                            tCommTable,
                                                                            tIdToIndMap,
                                                                            tMaxNumAdofs,
                                                                            tInterpolationMesh );

            if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
            {
                mModelSolverInterface->set_param("L2")= (sint)mDofOrder;
                mModelSolverInterface->set_param("TEMP")= (sint)mDofOrder;
            }

            //-----------------------------------------------------------------------------------------------------------

            moris::uint tEquationObjectCounter = tBlockSetCells.size();
            tElementBlockCounter = 0;
            mElementBlocks( tElementBlockCounter++ )->finalize( mModelSolverInterface );

            for( luint Ik = 0; Ik < aSidesetList.size(); ++Ik )
            {
                // get the treated sideset name
                std::string tTreatedSidesetName = tSideSetsNames( aSidesetList( Ik ) );

                // create a cell of sideset element
                moris::Cell< mtk::Cell const * > tSideSetElement;

                // create a list of sideset ordinal
                Matrix< IndexMat >  aSidesetOrdinals;

                // get the treated sideset elements and ordinals
                tInterpolationMesh->get_sideset_cells_and_ords( tTreatedSidesetName, tSideSetElement, aSidesetOrdinals );

                mElementBlocks( tElementBlockCounter++ )->finalize( mModelSolverInterface );

                for( luint k = 0; k < tSideSetElement.size(); ++k )
                {
                    mElements( tEquationObjectCounter )->set_list_of_side_ordinals( {{aSidesetOrdinals( k )}} ); //FIXME

                    // get the nodal weak bcs of the element
                    Matrix< DDRMat > & tNodalWeakBCs = mElements( tEquationObjectCounter )->get_weak_bcs();

                    // get the element number of nodes
                    uint tNumberOfNodes = mElements( tEquationObjectCounter++ )->get_num_nodes();

                    // set size of the element nodal weak bc
                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                    // loop over the element nodes
                    Matrix< IndexMat > tNodeIndices = tSideSetElement( k )->get_vertex_inds();

                    //--------------------------------------------------------------------------------------------
                    for( uint l = 0; l < tNumberOfNodes; l++ )
                    {
                        if ( aSidesetBCTypeList( Ik ) == fem::BC_Type::DIRICHLET )
                        {
                        // copy weak bc into element
                        tNodalWeakBCs( l ) = 5.0;
                        }
                        else if ( aSidesetBCTypeList( Ik ) == fem::BC_Type::NEUMANN )
                        {
                            // copy weak bc into element
                            tNodalWeakBCs( l ) = 20.0;
                        }
                    }
                }
            }

            //------------------------------------------------------------------------------------------
            mModelSolverInterface->finalize();

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created Model-Solver Interface in %5.3f seconds.\n\n",
                        ( double ) tElapsedTime / 1000 );
            }

        }

//------------------------------------------------------------------------------

        Model::~Model()
        {
            // delete SI
            delete mSolverInterface;

            // delete MSI
            delete mModelSolverInterface;

            // delete IWGs
            for( auto tIWG : mIWGs )
            {
                tIWG.clear();
            }

            // delete elements
            for( auto tElement : mElements )
            {
                delete tElement;
            }

            // delete nodes
            for( auto tNode : mNodes )
            {
                delete tNode;
            }
        }

//------------------------------------------------------------------------------

        void
        Model::set_weak_bcs( const Matrix<DDRMat> & aWeakBCs )
        {
            // set weak BCs
            for( auto tElement : mElements )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k < tNumberOfNodes; ++k )
                {
                    // copy weakbc into element
                    tNodalWeakBCs( k ) = aWeakBCs( tElement->get_node_index( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Model::set_weak_bcs_from_nodal_field( moris_index aFieldIndex )
        {
            for( auto tElement : mElements )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // copy weakbc into element
                    tNodalWeakBCs( k ) = mMeshManager->get_interpolation_mesh(0)->get_value_of_scalar_field( aFieldIndex,
                                                                           EntityRank::NODE,
                                                                           tElement->get_node_index( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        real
        Model::compute_integration_error(
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
        {
            real aError = 0.0;
            for( auto tElement : mElements )
            {
                aError += tElement->compute_integration_error( aFunction );
            }
            return aError;
        }

//------------------------------------------------------------------------------

        uint
        Model::get_lagrange_order_from_mesh()
        {

            // set order of this model according to Lagrange order
            // of first element on mesh
           return mtk::interpolation_order_to_uint(
            mMeshManager->get_interpolation_mesh(0)->get_mtk_cell( 0 ).get_interpolation_order() );
        }

//------------------------------------------------------------------------------

        void
        Model::set_dof_order( const uint aOrder )
        {
            mDofOrder = aOrder;
//            MORIS_ASSERT( aOrder == mDofOrder,
//                    "Model: the functionality to change the order of the model has nor been implemented yet" );
        }

//------------------------------------------------------------------------------
        real
        Model::compute_element_average( const uint aElementIndex )
        {
            return mElements( aElementIndex )->compute_element_average_of_scalar_field();
        }

//------------------------------------------------------------------------------

        void
        Model::output_solution( const std::string & aFilePath )
        {
            if ( mMeshManager->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR )
            {
                mSolHMR.set_size(mMeshManager->get_interpolation_mesh(0)->get_num_nodes(),1,-1.0);

                mtk::Interpolation_Mesh* tInterpMesh = mMeshManager->get_interpolation_mesh(0);

                moris::Cell<std::string> tBlockSetsNames = tInterpMesh->get_set_names( EntityRank::ELEMENT);

                for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
                {
                    Matrix< IndexMat > tBlockSetElementInd
                        = tInterpMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( Ik ) );

                    for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                    {
                       uint tNumVert = tInterpMesh->get_mtk_cell( k ).get_number_of_vertices();

                       //print( mElements(k)->get_pdof_values(), "Element");

                       for( luint Jk=0; Jk < tNumVert; ++Jk )
                       {
                           moris_index tID= tInterpMesh->get_mtk_cell( k ).get_vertex_pointers()( Jk) ->get_index();

                           mSolHMR(tID) = mElements(k)->get_pdof_values()(Jk);
                       }
                    }
                }
            }
            else
            {
            // 8) Postprocessing
            // dof type list for the solution to write on the mesh
            moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

            uint tNumOfNodes = mNodes.size();

            // create a matrix to be filled  with the solution
            Matrix< DDRMat > tTempSolutionField( tNumOfNodes, 1 );

            mtk::Interpolation_Mesh* tInterpMesh = mMeshManager->get_interpolation_mesh(0);

            // loop over the nodes
            for( uint i = 0; i < tNumOfNodes; i++ )
            {
                // get a list of elements connected to the ith node
                Matrix<IndexMat> tConnectedElements =
                        tInterpMesh->get_entity_connected_to_entity_loc_inds( static_cast< moris_index >( i ),
                                                                    EntityRank::NODE,
                                                                    EntityRank::ELEMENT );

                // number of connected element
                uint tNumConnectElem = tConnectedElements.numel();

                // reset the nodal value
                real tNodeVal = 0.0;

                // loop over the connected elements
                for( uint j = 0; j < tNumConnectElem; j++ )
                {
                    // extract the field value at the ith node for the jth connected element
                    real tElemVal = mElements( tConnectedElements( j ) )->get_element_nodal_pdof_value( i, tDofTypeList);
                    // add up the contribution of each element to the node value
                    tNodeVal = tNodeVal + tElemVal;
                }
                // fill the solution matrix with the node value
                tTempSolutionField( i ) = tNodeVal/tNumConnectElem;
            }
            //print( tTempSolutionField, "tTempSolutionField" );

            // add field to the mesh
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds( aFilePath,
                                                             EntityRank::NODE,
                                                             tTempSolutionField );

            // create output mesh
            std::string tOutputFile = "./int_ElemDiff_test_11.exo";
            tInterpMesh->create_output_mesh( tOutputFile );
            }
        }

        void
        Model::output_solution_nils_HACK( const std::string & aFilePath,
        		                                mtk::Mesh * aMesh )
        {
            if ( aMesh->get_mesh_type() == MeshType::HMR )
            {
                mSolHMR.set_size(aMesh->get_num_nodes(),1,-1.0);

                moris::Cell<std::string> tBlockSetsNames = aMesh->get_set_names( EntityRank::ELEMENT);

                for( luint Ik=0; Ik < tBlockSetsNames.size(); ++Ik )
                {
                    Matrix< IndexMat > tBlockSetElementInd
                        = aMesh->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetsNames( Ik ) );

                    for( luint k=0; k < tBlockSetElementInd.numel(); ++k )
                    {
                       uint tNumVert = aMesh->get_mtk_cell( k ).get_number_of_vertices();

                       //print( mElements(k)->get_pdof_values(), "Element");

                       for( luint Jk=0; Jk < tNumVert; ++Jk )
                       {
                           moris_index tID= aMesh->get_mtk_cell( k ).get_vertex_pointers()( Jk) ->get_index();

                           mSolHMR(tID) = mElements(k)->get_pdof_values()(Jk);
                       }
                    }
                }
            }
        }

    } /* namespace mdl */
} /* namespace moris */
