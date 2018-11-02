
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src


#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver_Manager.hpp"

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_print.hpp" // for check

// fixme: #ADOFORDERHACK
#include "MSI_Adof_Order_Hack.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

//        Model::Model(
//        		NLA::Nonlinear_Problem * aNonlinearProblem,
//        		std::shared_ptr< NLA::Nonlinear_Solver >  aSolver,
//                mtk::Mesh           * aMesh,
//                fem::IWG            & aIWG,
//                const Matrix< DDRMat > & aWeakBCs,
//                Matrix< DDRMat >       & aDOFs )

        Model::Model(
                mtk::Mesh           * aMesh,
                fem::IWG            * aIWG  ) : mMesh( aMesh )
        {

            // fixme: #ADOFORDERHACK
            if ( moris::MSI::gAdofOrderHack == 0 )
            {
                moris::MSI::gAdofOrderHack  = this->get_lagrange_order_from_mesh();
            }
            mDofOrder = moris::MSI::gAdofOrderHack;

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

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create elements
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            // ask mesh about number of elements on proc
            luint tNumberOfElements = aMesh->get_num_elems();

            // create equation objects
            mElements.resize( tNumberOfElements, nullptr );

            for( luint k=0; k<tNumberOfElements; ++k )
            {
                // create the element
                mElements( k ) = new fem::Element(
                        & aMesh->get_writable_mtk_cell( k ),
                        aIWG,
                        mNodes );
            }

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create Model Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mModelSolverInterface = new moris::MSI::Model_Solver_Interface(
                    mElements,
                    aMesh->get_communication_table(),
                    mCoefficientsMap,
                    aMesh->get_num_coeffs( mDofOrder ) );

            // calculate AdofMap
            mAdofMap = mModelSolverInterface->get_dof_manager()->get_adof_ind_map();

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: create Solver Interface
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            mSolverInterface =  new moris::MSI::MSI_Solver_Interface(
                    mModelSolverInterface,
                    mModelSolverInterface->get_dof_manager() );


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
            mNonlinerarSolver = tNonlinFactory.create_nonlinear_solver(
                    NLA::NonlinearSolverType::NEWTON_SOLVER );

            // create factory for linear solver
            dla::Solver_Factory  tSolFactory;

            // create linear solver
            mLinearSolver = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

            // set default parameters for linear solver
            mLinearSolver->set_param("AZ_diagnostics") = AZ_none;
            mLinearSolver->set_param("AZ_output") = AZ_none;

            // create solver manager
            mSolverManager = new dla::Linear_Solver_Manager();

            // set manager and settings
            mNonlinerarSolver->set_linear_solvers( mSolverManager );

            // set first solver
            mNonlinerarSolver->set_linear_solver( 0, mLinearSolver );
        }

//------------------------------------------------------------------------------

        Model::~Model()
        {

            // delete manager
            delete mSolverManager;

            // delete problem
            delete mNonlinerarProblem;

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
                    tNodalWeakBCs( k )
                            = mMesh->get_value_of_scalar_field(
                                    aFieldIndex,
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
            mNonlinerarSolver->solver_nonlinear_system( mNonlinerarProblem );

            // temporary array for solver
            Matrix< DDRMat > tSolution;
            mNonlinerarSolver->get_full_solution( tSolution );

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
            MORIS_ASSERT( aOrder == mDofOrder,
                    "Model: the functionality to change the order of the model has nor been implemented yet" );
        }
    } /* namespace mdl */
} /* namespace moris */
