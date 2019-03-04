
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "MTK_Tools.hpp"

#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"

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

//        Model::Model(
//                mtk::Mesh           * aMesh,
//                fem::IWG            * aIWG,
//                const uint    aBSplineOrder) : mMesh( aMesh )

        Model::Model(       mtk::Mesh * aMesh,
                      const uint        aBSplineOrder) : mMesh( aMesh )
        {

            // start timer
            tic tTimer1;

            mDofOrder = aBSplineOrder;

            if ( mDofOrder == 0 )
            {
                mDofOrder  = this->get_lagrange_order_from_mesh();
            }

            // get map from mesh
            mMesh->get_adof_map( mDofOrder, mCoefficientsMap );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create nodes
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // ask mesh about number of nodes on proc
            luint tNumberOfNodes = aMesh->get_num_nodes();

            // create node objects
            mNodes.resize(  tNumberOfNodes, nullptr );

            for( luint k=0; k<tNumberOfNodes; ++k )
            {
                mNodes( k ) = new fem::Node( &aMesh->get_mtk_vertex( k ) );
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

            //FIXME create IWGs- create L2 only
            // input a cell of IWG types to be created
            Cell< fem::IWG_Type > tIWGTypeList = { fem::IWG_Type::L2 };

            // number of IWGs to be created
            uint tNumOfIWGs = tIWGTypeList.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs for the problem considered
            Cell< fem::IWG* > tIWGs( tNumOfIWGs , nullptr );

            // loop over the IWG types
            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                // create an IWG with the factory for the ith IWG type
                tIWGs( i ) = tIWGFactory.create_IWGs( tIWGTypeList( i ) );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer2;

            // ask mesh about number of elements on proc
            luint tNumberOfElements = aMesh->get_num_elems();

            // create equation objects
            mElements.resize( tNumberOfElements, nullptr );

            for( luint k=0; k<tNumberOfElements; ++k )
            {
//                // create the element
//                mElements( k ) = new fem::Element( & aMesh->get_mtk_cell( k ),
//                                                   aIWG,
//                                                   mNodes );
                // create the element
                mElements( k ) = new fem::Element( & aMesh->get_mtk_cell( k ),
                                                     tIWGs,
                                                     mNodes );
            }

            if( par_rank() == 0)
            {
                // stop timer
                real tElapsedTime = tTimer2.toc<moris::chronos::milliseconds>().wall;

                // print output
                std::fprintf( stdout,"Model: created %u FEM elements in %5.3f seconds.\n\n",
                        ( unsigned int ) tNumberOfElements,
                        ( double ) tElapsedTime / 1000 );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create Model Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // start timer
            tic tTimer3;

            mModelSolverInterface = new moris::MSI::Model_Solver_Interface( mElements,
                                                                            aMesh->get_communication_table(),
                                                                            mCoefficientsMap,
                                                                            aMesh->get_num_coeffs( mDofOrder ),      //FIXME
                                                                            mMesh );

            mModelSolverInterface->set_param("L2")= (sint)mDofOrder;

            //mModelSolverInterface->set_param("L2") = (sint)aBSplineOrder;

            mModelSolverInterface->finalize();

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mSolverInterface =  new moris::MSI::MSI_Solver_Interface( mModelSolverInterface );


            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 5: create Nonlinear Problem
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mNonlinerarProblem =  new NLA::Nonlinear_Problem( mSolverInterface );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 6: create Solvers and solver manager
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // create factory for nonlinear solver
            NLA::Nonlinear_Solver_Factory tNonlinFactory;

            // create nonlinear solver
            mNonlinerarSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            // create factory for linear solver
            dla::Solver_Factory  tSolFactory;

            // create linear solver
            mLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

            // set default parameters for linear solver
            mLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            mLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
            //mLinearSolverAlgorithm->set_param("AZ_keep_info") = 1;
            //mLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;

            // create solver manager
            mLinSolver = new dla::Linear_Solver();
            mNonlinearSolver = new NLA::Nonlinear_Solver();

            // set manager and settings
            mNonlinerarSolverAlgorithm->set_linear_solver( mLinSolver );

            // set first solver
            mLinSolver->set_linear_algorithm( 0, mLinearSolverAlgorithm );

            mNonlinearSolver->set_nonlinear_algorithm( mNonlinerarSolverAlgorithm, 0 );

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

            // delete manager
            delete mLinSolver;

            // delete problem
            delete mNonlinerarProblem;

            // delete NonLinSolverManager
            delete mNonlinearSolver;

            // delete SI
            delete mSolverInterface;

            // delete MSI
            delete mModelSolverInterface;

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

                for( uint k=0; k<tNumberOfNodes; ++k )
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
                    tNodalWeakBCs( k ) = mMesh->get_value_of_scalar_field( aFieldIndex,
                                                                           EntityRank::NODE,
                                                                           tElement->get_node_index( k ) );
                }
            }
        }

//------------------------------------------------------------------------------

        void
        Model::solve( Matrix<DDRMat> & aSolution )
        {

            // call solver
            mNonlinearSolver->solve( mNonlinerarProblem );

            // temporary array for solver
            Matrix< DDRMat > tSolution;
            mNonlinerarSolverAlgorithm->get_full_solution( tSolution );

            // get length of array
            uint tLength = tSolution.length();

            // make sure that length of vector is correct
            MORIS_ERROR( tLength <= (uint)  mMesh->get_num_coeffs( mDofOrder ),
                    "Number of ADOFs does not match" );

            // rearrange data into output
            aSolution.set_size( tLength, 1 );

            for( uint k=0; k<tLength; ++k )
            {
                aSolution( k ) = tSolution( mAdofMap( k ) );
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
            mMesh->get_mtk_cell( 0 ).get_interpolation_order() );
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

    } /* namespace mdl */
} /* namespace moris */
