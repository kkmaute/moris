
// added by christian: link to Google Perftools
#ifdef WITHGPERFTOOLS
#include <gperftools/profiler.h>
#endif

#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src

#include "cl_Solver_Factory.hpp"
#include "cl_Solver_Input.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"

#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
//#include "cl_MSI_Node_Obj.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

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
            aMesh->get_adof_map( tAdofMap );

            // this part does not work yet in parallel
            auto tMSI = new moris::MSI::Model_Solver_Interface(
                    mElements,
                    aMesh->get_communication_table(),
                    tAdofMap,
                    aMesh->get_num_coeffs() );



            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Solver > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );


            // create interface
            moris::MSI::MSI_Solver_Interface *  tSolverInput;
            tSolverInput = new moris::MSI::MSI_Solver_Interface( tMSI, tMSI->get_dof_manager() );

            moris::Solver_Factory  tSolFactory;
            std::shared_ptr< Linear_Solver > tLin = tSolFactory.create_solver( tSolverInput, SolverType::AZTEC_IMPL );

            tNonLinSolver->set_linear_solver( tLin );

            tNonLinSolver->solver_nonlinear_system();


            Matrix< DDRMat > tDOFs;
            tNonLinSolver->get_full_solution( tDOFs );

            // -----------------
            uint tLength = tDOFs.length();

            // make sure that length of vector is correct
            MORIS_ASSERT( tLength == (uint)  aMesh->get_num_coeffs(),
                    "Number of ADOFs does not match" );

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

            // tidy up
            delete tSolverInput;

            // delete interface
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
