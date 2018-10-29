
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src

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
                fem::IWG            & aIWG,
                const Matrix< DDRMat > & aWeakBCs,
                Matrix< DDRMat >       & aDOFs )
        {
            // how many cells exist on current proc
            auto tNumberOfElements = aMesh->get_num_elems();

            // create nodes for these elements
            auto tNumberOfNodes = aMesh->get_num_nodes();

            // create node objects
            mNodes.resize(  tNumberOfNodes, nullptr );

            for( luint k=0; k<tNumberOfNodes; ++k )
            {
                mNodes( k ) = new fem::Node( &aMesh->get_mtk_vertex( k ) );
            }


            // create equation objects
            mElements.resize( tNumberOfElements, nullptr );

            for( luint k=0; k<tNumberOfElements; ++k )
            {
                // create the element
                mElements( k ) = new fem::Element(
                        & aMesh->get_writable_mtk_cell( k ), // <-- fixme does this need to be writable?
                        & aIWG,
                        mNodes,
                        aWeakBCs );
            }

            //for( luint k=0; k<tNumberOfElements; ++k )
           // {
           //     // compute matrix and RHS
           //     mElements( k )->compute_jacobian_and_residual();
           // }


            tic tTimer4;
            // create map for MSI
            map< moris_id, moris_index > tAdofMap;
            aMesh->get_adof_map( moris::MSI::gAdofOrderHack, tAdofMap ); // fixme: #ADOFORDERHACK

            // this part does not work yet in parallel
            auto tMSI = new moris::MSI::Model_Solver_Interface(
                    mElements,
                    aMesh->get_communication_table(),
                    tAdofMap,
                    aMesh->get_num_coeffs( moris::MSI::gAdofOrderHack ) ); // fixme: #ADOFORDERHACK

            moris::MSI::MSI_Solver_Interface *  tSolverInput;
            tSolverInput = new moris::MSI::MSI_Solver_Interface( tMSI, tMSI->get_dof_manager() );

            NLA::Nonlinear_Problem * tNonlinearProblem = new NLA::Nonlinear_Problem( tSolverInput );

            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Solver > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver > tLinSolver1 = tSolFactory.create_solver( SolverType::AZTEC_IMPL );


            tLinSolver1->set_param("AZ_diagnostics") = AZ_none;
            tLinSolver1->set_param("AZ_output") = AZ_none;

            dla::Linear_Solver_Manager * tLinearSolverManager = new  dla::Linear_Solver_Manager();

            //tLinearSolverManager->set_linear_solver( tLinSolver1 );
            tNonLinSolver->set_linear_solvers( tLinearSolverManager );
            tNonLinSolver->set_linear_solver( 0, tLinSolver1 );

            //tNonlinearProblem->set_interface( tSolverInput );

            tNonLinSolver->solver_nonlinear_system( tNonlinearProblem );
            Matrix< DDRMat > tDOFs;

            tNonLinSolver->get_full_solution( tDOFs );

            //print( tDOFs, "tdofs" );

            // -----------------
            uint tLength = tDOFs.length();

            // make sure that length of vector is correct
            MORIS_ERROR( tLength <= (uint)  aMesh->get_num_coeffs( moris::MSI::gAdofOrderHack ),
                    "Number of ADOFs does not match" );  // fixme: #ADOFORDERHACK

            // fixme this is only temporary. Needed for integration error
            //for( auto tElement : mElements )
            //{
            //     tElement->extract_values( tLin );
            //}

            auto tMap = tMSI->get_dof_manager()->get_adof_ind_map();

            aDOFs.set_size( tLength, 1 );
            for( uint k=0; k<tLength; ++k )
            {
                //aDOFs( tMap( k ) ) = tDOFs( k );
                aDOFs( k ) = tDOFs( tMap( k ) );
            }

            delete tLinearSolverManager;

            // tidy up
            delete tSolverInput;

            //delete interface
            delete tMSI;



        }

//------------------------------------------------------------------------------

        Model::~Model()
        {

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

    } /* namespace mdl */
} /* namespace moris */
