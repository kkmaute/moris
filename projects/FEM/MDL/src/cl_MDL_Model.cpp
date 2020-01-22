
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"              //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"       //MTK/src

#include "cl_FEM_Enums.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Model.hpp"

//#include "cl_DLA_Solver_Factory.hpp"
//#include "cl_DLA_Solver_Interface.hpp"

//#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
//#include "cl_NLA_Nonlinear_Solver.hpp"
//#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
//#include "cl_DLA_Linear_Solver_Aztec.hpp"
//#include "cl_DLA_Linear_Solver.hpp"

//#include "cl_TSA_Time_Solver_Factory.hpp"
//#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_VIS_Output_Manager.hpp"

#include "op_equal_equal.hpp"

#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------
    Model::Model(       mtk::Mesh_Manager                 * aMeshManager,
                  const uint                                aBSplineIndex,
                        moris::Cell< fem::Set_User_Info > & aSetInfo,
                  const moris_index                         aMeshPairIndex,
                  const bool                                aUseMultigrid ) : mMeshManager( aMeshManager ),
                                                                              mMeshPairIndex( aMeshPairIndex ),
                                                                              mUseMultigrid( aUseMultigrid )
    {
        mFemModel =  std::make_shared< fem::FEM_Model >( aMeshManager, aMeshPairIndex, aSetInfo );

        // FIXME comment???
        mBSplineIndex = aBSplineIndex;

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 0: initialize
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Get pointers to interpolation and integration mesh
        moris::moris_index tMeshPairIndex = aMeshPairIndex;
        mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
        mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
        mMeshManager->get_mesh_pair( tMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Create Model Solver Interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        // start timer
        tic tTimer1;

        //--------------------------FIXME------------------------------------
        // This part should not be needed anymore when MTK has all the functionalities
        Matrix< IdMat > tCommTable;
        moris::map< moris::moris_id, moris::moris_index > tIdToIndMap;
        moris::uint tMaxNumAdofs;

        if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
        {
            // get map from mesh
            tInterpolationMesh->get_adof_map( mBSplineIndex, mCoefficientsMap );

            tCommTable   = tInterpolationMesh->get_communication_table();
            tIdToIndMap  = mCoefficientsMap;
            tMaxNumAdofs = tInterpolationMesh->get_num_coeffs( mBSplineIndex );
        }
        else
        {
            tCommTable.set_size( 1, 1, 0 );
            tMaxNumAdofs = 15000000;
        }
        //--------------------------END FIXME--------------------------------

        mEquationSets = mFemModel->get_equation_sets();

        mModelSolverInterface = new moris::MSI::Model_Solver_Interface( mEquationSets,
                                                                        tCommTable,
                                                                        tIdToIndMap,
                                                                        tMaxNumAdofs,
                                                                        tInterpolationMesh );

        if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
        {
            mModelSolverInterface->set_param("L2")= (sint)mBSplineIndex;
            mModelSolverInterface->set_param("TEMP")= (sint)mBSplineIndex;
        }

        //------------------------------------------------------------------------------

        // finalize the fem sets
        for( luint Ik = 0; Ik < mEquationSets.size(); ++Ik )
        {
            // finalize the fem sets
        	mEquationSets( Ik )->finalize( mModelSolverInterface );
        }

        mModelSolverInterface->finalize( mUseMultigrid );

        // calculate AdofMap
        mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: create Solver Interface
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

        mSolverInterface->set_model( this );

        if( par_rank() == 0)
        {
            // stop timer
            real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

            // print output
            MORIS_LOG_INFO( "Model: created Model-Solver Interface in %5.3f seconds.\n\n",
                    ( double ) tElapsedTime / 1000 );
        }
    }

//------------------------------------------------------------------------------

    Model::Model(       mtk::Mesh_Manager*                  aMeshManager,
                  const uint                                aBSplineIndex,
                        moris::Cell< fem::Set_User_Info > & aSetInfo,
                        MSI::Design_Variable_Interface * aDesignVariableInterface,
                  const moris_index                         aMeshPairIndex,
                  const bool                                aUseMultigrid ) : mMeshManager( aMeshManager ),
                                                                              mMeshPairIndex( aMeshPairIndex ),
                                                                              mUseMultigrid( aUseMultigrid )
        {
            mFemModel =  std::make_shared< fem::FEM_Model >( aMeshManager, aMeshPairIndex, aSetInfo );

            // start timer
            tic tTimer1;

            // FIXME comment???
            mBSplineIndex = aBSplineIndex;

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 0: initialize
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Get pointers to interpolation and integration mesh
            moris::moris_index tMeshPairIndex = aMeshPairIndex;
            mtk::Interpolation_Mesh* tInterpolationMesh = nullptr;
            mtk::Integration_Mesh*   tIntegrationMesh   = nullptr;
            mMeshManager->get_mesh_pair( tMeshPairIndex, tInterpolationMesh, tIntegrationMesh );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // Create Model Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            //--------------------------FIXME------------------------------------
            // This part should not be needed anymore when MTK has all the functionalities
            Matrix< IdMat > tCommTable;
            moris::map< moris::moris_id, moris::moris_index > tIdToIndMap;
            moris::uint tMaxNumAdofs;

            if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
            {
                // get map from mesh
                tInterpolationMesh->get_adof_map( mBSplineIndex, mCoefficientsMap );

                tCommTable   = tInterpolationMesh->get_communication_table();
                tIdToIndMap  = mCoefficientsMap;
                tMaxNumAdofs = tInterpolationMesh->get_num_coeffs( mBSplineIndex );
            }
            else
            {
                tCommTable.set_size( 1, 1, 0 );
                tMaxNumAdofs = 15000000;
            }
            //--------------------------END FIXME--------------------------------

            mEquationSets = mFemModel->get_equation_sets();

            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( mEquationSets,
                                                                            tCommTable,
                                                                            tIdToIndMap,
                                                                            tMaxNumAdofs,
                                                                            tInterpolationMesh );

            if ( tInterpolationMesh->get_mesh_type() == MeshType::HMR )
            {
                mModelSolverInterface->set_param("L2")= (sint)mBSplineIndex;
                mModelSolverInterface->set_param("TEMP")= (sint)mBSplineIndex;
            }

            //------------------------------------------------------------------------------

            // finalize the fem sets
            for( luint Ik = 0; Ik < mEquationSets.size(); ++Ik )
            {
            	mEquationSets( Ik )->set_Dv_interface( aDesignVariableInterface);
                // finalize the fem sets
            	mEquationSets( Ik )->finalize( mModelSolverInterface );
            }

            mModelSolverInterface->finalize( mUseMultigrid );

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );

            mSolverInterface->set_model( this );

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer1.toc<moris::chronos::milliseconds>().wall;

                // print output
                MORIS_LOG_INFO( "Model: created Model-Solver Interface in %5.3f seconds.\n\n",
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

            mEquationSets.clear();

            // delete the fem cluster
            mEquationObjects.clear();
        }

//------------------------------------------------------------------------------

        map< moris_index, moris_index > & Model::get_mesh_set_to_fem_set_index_map( )
        {
            return mFemModel->get_mesh_set_to_fem_set_index_map();
        }

//------------------------------------------------------------------------------

        void Model::set_weak_bcs( const Matrix<DDRMat> & aWeakBCs )
        {
        	moris::Cell< MSI::Equation_Object * > tFemClusters = mFemModel->get_equation_objects();
            // set weak BCs
            for( auto tElement : tFemClusters )
            {
                Matrix< DDRMat > & tNodalWeakBCs = tElement->get_weak_bcs();
                uint tNumberOfNodes = tElement->get_num_nodes();
                tNodalWeakBCs.set_size( tNumberOfNodes, 1 );

                for( uint k=0; k < tNumberOfNodes; ++k )
                {
                    // copy weak bc into element
                    tNodalWeakBCs( k ) = aWeakBCs( tElement->get_node_index( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void Model::set_weak_bcs_from_nodal_field( moris_index aFieldIndex )
        {
        	moris::Cell< MSI::Equation_Object * > tFemClusters = mFemModel->get_equation_objects();

            for( auto tElement : tFemClusters )
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

        real Model::compute_integration_error(
                real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
        {
            real aError = 0.0;
            for( auto tElement : mEquationObjects )
            {
                aError += tElement->compute_integration_error( aFunction );
            }
            return aError;
        }

//------------------------------------------------------------------------------

        uint Model::get_lagrange_order_from_mesh()
        {

            // set order of this model according to Lagrange order
            // of first element on mesh
           return mtk::interpolation_order_to_uint( mMeshManager->get_interpolation_mesh(0)
                                                                ->get_mtk_cell( 0 ).get_interpolation_order() );
        }

//------------------------------------------------------------------------------

        void Model::set_dof_order( const uint aBSplineIndex )
        {
            mBSplineIndex = aBSplineIndex;
//            MORIS_ASSERT( aOrder == mBSplineIndex,
//                    "Model: the functionality to change the order of the model has nor been implemented yet" );
        }

//------------------------------------------------------------------------------
        real Model::compute_element_average( const uint aElementIndex )
        {
            return mEquationObjects( aElementIndex )->compute_element_average_of_scalar_field();
        }

//------------------------------------------------------------------------------

        void Model::output_solution( const std::string & aFilePath )
        {
            mSolHMR.set_size(mMeshManager->get_interpolation_mesh(0)->get_num_nodes(),1,0.0);

            mtk::Interpolation_Mesh* tInterpMesh = mMeshManager->get_interpolation_mesh(0);

            moris::Cell<std::string> tBlockSetsNames = tInterpMesh->get_set_names( EntityRank::ELEMENT);

            // iterate through fem blocks
            uint tNumFemSets = mEquationSets.size();

            for(moris::uint iSet = 0; iSet<tNumFemSets; iSet++)
            {
                // access set
                MSI::Equation_Set * tFemSet = mEquationSets(iSet);

                // access the element type in the set
                enum fem::Element_Type tElementTypeInSet = reinterpret_cast< fem::Set * >(tFemSet)->get_element_type();

                // if we are a bulk element, then we output the solution
                if(tElementTypeInSet == fem::Element_Type::BULK)
                {
                    // get mtk cell clusters in set
                    moris::Cell< mtk::Cluster const* > tCellClustersInSet = reinterpret_cast< fem::Set * >(tFemSet)->get_clusters_on_set();

                    // get coressponding equation objects
                    moris::Cell< MSI::Equation_Object * > & tEquationObj = reinterpret_cast< fem::Set * >(tFemSet)->get_equation_object_list();

                    // iterate through clusters in set
                    for(moris::uint iCl = 0; iCl < tCellClustersInSet.size(); iCl++)
                    {
                        mtk::Cluster const * tCluster = tCellClustersInSet(iCl);

                        moris::mtk::Cell const & tInterpCell = tCluster->get_interpolation_cell();

                        uint tNumVert = tInterpCell.get_number_of_vertices();

                        Matrix< DDRMat > & tPDofVals = tEquationObj(iCl)->get_pdof_values();

                        moris::Cell<moris::mtk::Vertex *> tVertices = tInterpCell.get_vertex_pointers();

                        for( luint Jk=0; Jk < tNumVert; ++Jk )
                        {
                            mSolHMR(tVertices(Jk)->get_index()) = tPDofVals(Jk);
                        }
                    }
                }

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

//            if ( mMeshManager->get_interpolation_mesh(0)->get_mesh_type() == MeshType::HMR )
//            {

                // iterate through fem blocks
                uint tNumFemSets = mEquationSets.size();

                for(moris::uint iSet = 0; iSet<tNumFemSets; iSet++)
                {
                    // access set
                    MSI::Equation_Set * tFemSet = mEquationSets(iSet);

                    // access the element type in the set
                    enum fem::Element_Type tElementTypeInSet = reinterpret_cast< fem::Set * >(tFemSet)->get_element_type();

                    // if we are a bulk element, then we output the solution
                    if(tElementTypeInSet == fem::Element_Type::BULK)
                    {
                        // access dof list
                        moris::Cell< enum MSI::Dof_Type > & tUniqueDofList = reinterpret_cast< fem::Set * >(tFemSet)->get_unique_dof_type_list();

                        // see if the dof of interest is in there
                        bool tDofInSet = dof_type_is_in_list(aDofType,tUniqueDofList);

                        // if the dof is in this set, we need to add solution information to output
                        if(tDofInSet)
                        {
                            // get mtk cell clusters in set
                            moris::Cell< mtk::Cluster const* > tCellClustersInSet = reinterpret_cast< fem::Set * >(tFemSet)->get_clusters_on_set();

                            // get coressponding equation objects
                            moris::Cell< MSI::Equation_Object * > & tEquationObj = reinterpret_cast< fem::Set * >(tFemSet)->get_equation_object_list();

                            enum MSI::Dof_Type tDofType = MSI::Dof_Type::UX;
                            if ( aDofType == MSI::Dof_Type::TEMP)
                            {
                                tDofType = MSI::Dof_Type::TEMP;
                            }

                            // get the field interpolator
                            //fem::Field_Interpolator_Manager * tFIM = tFemSet->get_field_interpolators_manager();
                            fem::Field_Interpolator_Manager * tFIM = reinterpret_cast< fem::Set * >(tFemSet)->get_field_interpolator_manager( mtk::Master_Slave::MASTER );
                            fem::Field_Interpolator* tFieldInterp = tFIM->get_field_interpolators_for_type(tDofType);

                            // iterate through clusters in set
                            for(moris::uint iCl = 0; iCl < tCellClustersInSet.size(); iCl++)
                            {
                                // get the cluster
                                mtk::Cluster const * tCluster = tCellClustersInSet(iCl);

                                tEquationObj(iCl)->compute_my_pdof_values( );

                                // get the pdof values for this cluster
                                Cell < Matrix< DDRMat > >  tPDofVals;
                                tEquationObj(iCl)->get_my_pdof_values( tFieldInterp->get_dof_type(), tPDofVals);

                                //---------------------------------------------------------------------------------
                                uint tCols = tPDofVals.size();
                                uint tRows = tPDofVals( 0 ).numel();

                                Matrix< DDRMat > tReshapedPdofValues( tRows, tCols );

                                for( uint Ik = 0; Ik < tCols; Ik++ )
                                {
                                    tReshapedPdofValues( { 0, tRows - 1 }, { Ik, Ik } ) = tPDofVals( Ik ).matrix_data();
                                }

                                moris_index tDofIndex = 0;

                                for( uint Ik = 0; Ik < tFieldInterp->get_dof_type().size(); Ik ++ )
                                {
                                    if( aDofType == tFieldInterp->get_dof_type()( Ik ) )
                                    {
                                        tDofIndex = Ik;
                                    }
                                }
                                //---------------------------------------------------------------------------------

                                // set coefficients in field interpolator
                                tFieldInterp->set_coeff(tReshapedPdofValues);



                                // check if its trivial
                                bool tTrivialCluster = tCluster->is_trivial();

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
                                         tSolutionOnInteg(tVertIndex) = tPDofVals( tDofIndex )(iVert);

                                         // update count that this vertex has been encountered
                                         tVertexCount(tVertIndex) = tVertexCount(tVertIndex) + 1.0;
                                    }

                                }
                                else
                                {
                                    // get the vertices in the cluster
                                    moris::Cell<moris::mtk::Vertex *> tVertices = tCluster->get_primary_vertices_in_cluster();

                                    // iterate through vertices in cluster
                                    for(moris::uint iVert =0; iVert<tVertices.size(); iVert++)
                                    {
                                        // get vert spatial parametric coords
                                        moris::Matrix<moris::DDRMat> tVertParamCoords = tCluster->get_vertex_local_coordinate_wrt_interp_cell(tVertices(iVert));

                                        // fixme: figure out how to not need this transpose
                                        if(isrow(tVertParamCoords))
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
                                        tSolutionOnInteg(tVertIndex) = tSolFieldAtIntegPoint(tDofIndex);

                                        // update count that this vertex has been encountered
                                        tVertexCount(tVertIndex) = tVertexCount(tVertIndex) + 1.0;
                                    }

                                }

                            }

                        }
                    }
                }

                // iterate through and divide by the number of times a vertex has been added to solution field
//                for(moris::uint iVert =0; iVert<tNumVertsInIGMesh; iVert++)
//                {
//                    if(tVertexCount(iVert) > 0.0)
//                    {
//                        tSolutionOnInteg(iVert) = tSolutionOnInteg(iVert)/tVertexCount(iVert);
//                    }
//                }

//            }

            return tSolutionOnInteg;
        }

        bool Model::dof_type_is_in_list( enum MSI::Dof_Type aDofTypeToFind,
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

//------------------------------------------------------------------------------

        void Model::output_solution( const uint aVisMeshIndex,
                                     const real aTime )
        {
            // create the visualization mesh for this index
            mOutputManager->create_visualization_mesh( aVisMeshIndex,
                                                       mMeshManager,
                                                       mMeshPairIndex );

            // set the visualization sets to the equation sets
            mOutputManager->set_visualization_sets( aVisMeshIndex,
                                                    this );

            // Write the mesh to file
            mOutputManager->write_mesh( aVisMeshIndex,
                                        aTime );

            // write requested fields
            mOutputManager->write_field( aVisMeshIndex,
                                         this );

            // end writing and delete vis mesh
            mOutputManager->end_writing( aVisMeshIndex );
        }

    } /* namespace mdl */
} /* namespace moris */
