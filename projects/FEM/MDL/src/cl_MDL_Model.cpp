
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"              //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"       //MTK/src

#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element_Bulk.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_Set.hpp"

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

#include "op_equal_equal.hpp"

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_print.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"


namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        Model::Model(       mtk::Mesh_Manager*                          aMeshManager,
                      const uint                                        aBSplineOrder,
                            moris::Cell< moris::Cell< fem::IWG_Type > > aIWGTypeList,
                      const moris::Cell< moris_index >                  & aBlocksetList,
                      const moris::Cell< moris_index >                  aSidesetList,
                      const moris::Cell< fem::BC_Type >                 aSidesetBCTypeList,
                      const moris::Cell< moris_index >                  aDoubleSidesetList,
                      const moris_index                                 aMeshPairIndex) : mMeshManager( aMeshManager ),
                                                                                          mMeshPairIndex( aMeshPairIndex )
        {
            // start timer
            tic tTimer1;

            mDofOrder = aBSplineOrder;


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: initialize
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Get pointers to interpolation and integration mesh
            moris::moris_index tMeshPairIndex = aMeshPairIndex;
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMeshManager->get_mesh_pair( tMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

//            MORIS_ERROR( !(mDofOrder == 0 && tInterpolationMesh->get_mesh_type() == MeshType::HMR), " HMR B-Spline order can't be 0");

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // ask mesh about number of nodes on proc
            luint tNumOfIPNodes = tInterpolationMesh->get_num_nodes();

            // create node objects
            mIPNodes.resize(  tNumOfIPNodes, nullptr );

            for( luint k = 0; k < tNumOfIPNodes; ++k )
            {
                mIPNodes( k ) = new fem::Node( &tInterpolationMesh->get_mtk_vertex( k ) );
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout, "Model: created %u FEM IP nodes in %5.3f seconds.\n\n",
                        ( unsigned int ) tNumOfIPNodes,
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

            // init the number of set
            uint tNumberOfFemSet = aBlocksetList.size() + aSidesetList.size() + aDoubleSidesetList.size();

            // create equation objects
            mFemSets.resize( tNumberOfFemSet, nullptr );
            mFemClusters.reserve( tNumberOfEquationObjects );

            //  Create Blockset Elements ---------------------------------------------------
            std::cout<<" Create Block-set Elements "<<std::endl;
            //------------------------------------------------------------------------------
            // init the fem set counter
            moris::uint tFemSetCounter = 0;

            // loop over the used mesh block-set
            for( luint Ik = 0; Ik < aBlocksetList.size(); ++Ik )
            {
                // create a list of cell clusters (this needs to stay in scope somehow)
                moris::mtk::Set * tBlockSet = tIntegrationMesh->get_block_by_index( aBlocksetList( Ik ) );
//                moris::Cell<mtk::Cluster const*> tBlockSetClusterList = tIntegrationMesh->get_block_by_index( aBlocksetList( Ik ) )
//                                                                                        ->get_clusters_on_set();

                // create new fem set
                mFemSets( tFemSetCounter ) = new fem::Set( tBlockSet,
                                                           fem::Element_Type::BULK,
                                                           mIWGs( 0 ),
                                                           mIPNodes );

                // collect equation objects associated with the block-set
                mFemClusters.append( mFemSets( tFemSetCounter )->get_equation_object_list() );

                // update fem set counter
                tFemSetCounter++;
            }

            // FIXME for setting the BC in the model, will be ultimately removed
            moris::uint tEquationObjectCounter = mFemClusters.size();

            //  Create Sideset Elements ---------------------------------------------------
            std::cout<<" Create Side-set Elements "<<std::endl;
            //------------------------------------------------------------------------------
            // loop over the used mesh side-set
            for( luint Ik = 0; Ik < aSidesetList.size(); ++Ik )
            {
                // create a list of side clusters
                moris::mtk::Set * tSideSet = tIntegrationMesh->get_side_set_by_index( aSidesetList( Ik ) );
//                moris::Cell<mtk::Cluster const*> tSideSetClusterList = tIntegrationMesh->get_side_set_by_index( aSidesetList( Ik ) )
//                                                                                        ->get_clusters_on_set();

                // create a new fem set
                mFemSets( tFemSetCounter ) = new fem::Set( tSideSet,
                                                           fem::Element_Type::SIDESET,
                                                           mIWGs( Ik + 1 ), //FIXME this is why we had a problem with Dirichlet and Neumann order
                                                           mIPNodes );

                // collect equation objects associated with the side-set
                mFemClusters.append( mFemSets( tFemSetCounter )->get_equation_object_list() );

                // update fem set counter
                tFemSetCounter++;
            }

            //  Create Double Side-set Elements ---------------------------------------------
            std::cout<<" Create Double Side-set Elements "<<std::endl;
            //------------------------------------------------------------------------------

            // loop over the used mesh double side-set
            for( luint Ik = 0; Ik < aDoubleSidesetList.size(); ++Ik )
            {
                // create a list of double side clusters
                moris::mtk::Set * tDoubleSideSet = tIntegrationMesh->get_double_side_set_by_index( aDoubleSidesetList( Ik ) );
//                moris::Cell<mtk::Cluster const*> tDoubleSideSetClusterList = tIntegrationMesh->get_double_side_set_by_index( aDoubleSidesetList( Ik ) )
//                                                                                        ->get_clusters_on_set();

                // create a new fem set
                mFemSets( tFemSetCounter ) = new fem::Set( tDoubleSideSet,
                                                           fem::Element_Type::DOUBLE_SIDESET,
                                                           mIWGs( 3 ), //FIXME this is why we had a problem with Dirichlet and Neumann order
                                                           mIPNodes );

                // collect equation objects associated with the side-set
                mFemClusters.append( mFemSets( tFemSetCounter )->get_equation_object_list() );

                // update fem set counter
                tFemSetCounter++;
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output

                std::fprintf( stdout,"Model: created %u FEM elements in %5.3f seconds.\n\n",
                        ( unsigned int ) mFemClusters.size(),
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
            }
            else
            {
                tCommTable.set_size( 1, 1, 0 );
                tMaxNumAdofs = 1000000;
            }
            //--------------------------END FIXME--------------------------------

            // Construct cell of equation sets (cast from child to base class of entire cell)
            moris::Cell< MSI::Equation_Set * > tEquationSets(mFemSets.size());
            for(moris::uint iSet = 0; iSet < mFemSets.size(); iSet++ )
            {
                tEquationSets(iSet) = mFemSets(iSet);
            }


            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( tEquationSets,
                                                                            tCommTable,
                                                                            tIdToIndMap,
                                                                            tMaxNumAdofs,
                                                                            tInterpolationMesh );

            if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
            {
                mModelSolverInterface->set_param("L2")= (sint)mDofOrder;
                mModelSolverInterface->set_param("TEMP")= (sint)mDofOrder;
            }

            //------------------------------------------------------------------------------

            // finalize the fem sets
            for( luint Ik = 0; Ik < mFemSets.size(); ++Ik )
            {
                // finalize the fem sets
                mFemSets( Ik )->finalize( mModelSolverInterface );
            }

            //------------------------------------------------------------------------------
            // FIXME Imposing boundary conditions
            //------------------------------------------------------------------------------

            // get the side-set names
            moris::Cell<std::string> tSidesetNames = tIntegrationMesh->get_set_names( EntityRank::FACE );

            // loop over the side-set cluster
            for( luint Ik = 0; Ik < aSidesetList.size(); ++Ik )
            {
                // get the treated sideset name
                std::string tTreatedSidesetName = tSidesetNames( aSidesetList( Ik ) );

                // get the sideset ordinal
                moris_index tSideSetOrd = tIntegrationMesh->get_side_set_index( tTreatedSidesetName );

                // get the side-sets clusters
                moris::Cell<moris::mtk::Cluster const *> tClustersInSideSet = tIntegrationMesh->get_side_set_cluster(tSideSetOrd);

                for( luint k = 0; k < tClustersInSideSet.size(); ++k )
                {
                    // get the nodal weak bcs of the element
                    Matrix< DDRMat > & tNodalWeakBCs = mFemClusters( tEquationObjectCounter )->get_weak_bcs();

                    // get the element number of nodes
                    uint tNumberOfNodes = mFemClusters( tEquationObjectCounter++ )->get_num_nodes();

                    // set size of the element nodal weak bc
                    tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

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
            mIWGs.clear();

            // delete fem nodes
            mIPNodes.clear();

            // delete the fem sets
            mFemSets.clear();

            // delete the fem cluster
            mFemClusters.clear();
        }

//------------------------------------------------------------------------------

        void
        Model::set_weak_bcs( const Matrix<DDRMat> & aWeakBCs )
        {
            // set weak BCs
            for( auto tElement : mFemClusters )
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
            for( auto tElement : mFemClusters )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k<tNumberOfNodes; ++k )
                {
                    // copy weakbc into element
                    tNodalWeakBCs( k ) = mMeshManager->get_interpolation_mesh( mMeshPairIndex )
                                                     ->get_value_of_scalar_field( aFieldIndex,
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
            for( auto tElement : mFemClusters )
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
            return mFemClusters( aElementIndex )->compute_element_average_of_scalar_field();
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

                           mSolHMR(tID) = mFemClusters(k)->get_pdof_values()(Jk);
                       }
                    }
                }
            }
            else
            {
            // 8) Postprocessing
            // dof type list for the solution to write on the mesh
            moris::Cell< MSI::Dof_Type > tDofTypeList = { MSI::Dof_Type::TEMP };

            uint tNumOfNodes = mIPNodes.size();

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
                    real tElemVal = mFemClusters( tConnectedElements( j ) )->get_element_nodal_pdof_value( i, tDofTypeList);
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

        Matrix<DDRMat>
        Model::get_solution_for_integration_mesh_output( enum MSI::Dof_Type aDofType )
        {
            // number of vertices in integration mesh
            uint tNumVertsInIGMesh = mMeshManager->get_integration_mesh(mMeshPairIndex)->get_num_entities(EntityRank::NODE);

            //  initialize integration mesh solution
            Matrix<DDRMat> tSolutionOnInteg ( tNumVertsInIGMesh,1,0.0);

            // keep track of how many times I have added to a given node
            Matrix<DDRMat> tVertexCount ( tNumVertsInIGMesh,1,0.0);

            if ( mMeshManager->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR )
            {

                // iterate through fem blocks
                uint tNumFemSets = mFemSets.size();

                for(moris::uint iSet = 0; iSet<tNumFemSets; iSet++)
                {
                    // access set
                    fem::Set * tFemSet = mFemSets(iSet);

                    // access the element type in the set
                    enum fem::Element_Type tElementTypeInSet = tFemSet->get_set_element_type();

                    // if we are a bulk element, then we output the solution
                    if(tElementTypeInSet == fem::Element_Type::BULK)
                    {
                        // access dof list
                        moris::Cell< enum MSI::Dof_Type > & tUniqueDofList = tFemSet->get_unique_dof_type_list();

                        // see if the dof of interest is in there
                        bool tDofInSet = dof_type_is_in_list(aDofType,tUniqueDofList);

                        // if the dof is in this set, we need to add solution information to output
                        if(tDofInSet)
                        {
                            // get mtk cell clusters in set
                            moris::Cell< mtk::Cluster const* > tCellClustersInSet = tFemSet->get_clusters_on_set();

                            // get coressponding equation objects
                            moris::Cell< MSI::Equation_Object * > & tEquationObj = tFemSet->get_equation_object_list();

                            // get the field interpolator
                            fem::Field_Interpolator* tFieldInterp = tFemSet->get_dof_type_field_interpolators(aDofType);

                            // iterate through clusters in set
                            for(moris::uint iCl = 0; iCl < tCellClustersInSet.size(); iCl++)
                            {
                                // get the cluster
                                mtk::Cluster const * tCluster = tCellClustersInSet(iCl);

                                // get the pdof values for this cluster
                                Matrix< DDRMat > & tPDofVals = tEquationObj(iCl)->get_pdof_values();

                                tFieldInterp->set_coeff(tPDofVals);

                                // check if its trivial
                                bool tTrivialCluster = tCluster->is_trivial();

                                // get the vertices in the cluster
                                moris::Cell<moris::mtk::Vertex const *> const & tVertices = tCluster->get_vertices_in_cluster();

                                // in the trivial case, the interp pdofs, correspond to the integration vertices in clusters
                                if(tTrivialCluster)
                                {
                                    // get integration cell
                                    moris::Cell<moris::mtk::Cell const *> const & tPrimaryCells = tCluster->get_primary_cells_in_cluster();

                                    MORIS_ASSERT(tPrimaryCells.size() , "There needs to be exactly 1 primary cell in cluster for the trivial case");

                                    // an assumption is made here that the vertices on the primary correspond to the index in pdof vector
                                    moris::Cell< moris::mtk::Vertex* > tVerticesOnPrimaryCell = tPrimaryCells(0)->get_vertex_pointers();

                                    // iterate through vertices and add their pdof value to integration mesh solution vector
                                    for(moris::uint iVert =0; iVert<tVerticesOnPrimaryCell.size(); iVert++)
                                    {
                                        // vertex index
                                         moris_index tVertIndex = tVerticesOnPrimaryCell(iVert)->get_index();

                                         // add to solution vector
                                         tSolutionOnInteg(tVertIndex) = tSolutionOnInteg(tVertIndex) + tPDofVals(iVert);

                                         // update count that this vertex has been encountered
                                         tVertexCount(tVertIndex) = tVertexCount(tVertIndex) + 1.0;
                                    }

                                }
                                else
                                {
                                    // iterate through vertices in cluster
                                    for(moris::uint iVert =0; iVert<tVertices.size(); iVert++)
                                    {
                                        // get vert spatial parametric coords
                                        moris::Matrix<moris::DDRMat> tVertParamCoords = tCluster->get_vertex_local_coordinate_wrt_interp_cell(tVertices(iVert));

                                        // fixme: figure out how to not need this transpose
                                        if(iscol(tVertParamCoords))
                                        {
                                            tVertParamCoords = trans(tVertParamCoords);
                                        }

                                        // size
                                        uint tSpatialParamSize = tVertParamCoords.numel();

                                        // add a time parametric coord
                                        tVertParamCoords.resize(tSpatialParamSize+1,1);
                                        tVertParamCoords(tSpatialParamSize) = 0;

                                        // set in field interpolator
                                        tFieldInterp->set_space_time(tVertParamCoords);

                                        // evaluate field value at this point
                                        Matrix<DDRMat> tSolFieldAtIntegPoint = tFieldInterp->val();

                                        // vertex index
                                        moris_index tVertIndex = tVertices(iVert)->get_index();

                                        // add to solution vector
                                        tSolutionOnInteg(tVertIndex) = tSolutionOnInteg(tVertIndex) + tSolFieldAtIntegPoint(0);

                                        // update count that this vertex has been encountered
                                        tVertexCount(tVertIndex) = tVertexCount(tVertIndex) + 1.0;
                                    }

                                }

                            }

                        }
                    }
                }

                // iterate through and divide by the number of times a vertex has been added to solution field
                for(moris::uint iVert =0; iVert<tNumVertsInIGMesh; iVert++)
                {
                    if(tVertexCount(iVert) > 0.0)
                    {
                        tSolutionOnInteg(iVert) = tSolutionOnInteg(iVert)/tVertexCount(iVert);
                    }
                }

            }

            return tSolutionOnInteg;
        }

        bool
        Model::dof_type_is_in_list( enum MSI::Dof_Type aDofTypeToFind,
                                    moris::Cell< enum MSI::Dof_Type > & aDofList)
        {
            bool tDofInList = false;

            uint tNumDofsInFullList = aDofList.size();

            // iterate through
            for(moris::uint iDof = 0; iDof<tNumDofsInFullList; iDof++)
            {
                if(aDofList(iDof) == aDofTypeToFind)
                {
                    tDofInList  = true;
                    break;
                }
            }
            return tDofInList;
        }

    } /* namespace mdl */
} /* namespace moris */
