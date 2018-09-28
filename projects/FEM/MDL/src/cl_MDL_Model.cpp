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

#include "cl_MSI_Model_Solver_Interface.hpp"

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        Model::Model(
                mtk::Mesh              * aMesh,
                fem::IWG               & aIWG,
                const Matrix< DDRMat > & aWeakBCs,
                Matrix< DDRMat >       & aDOFs )
        {
            // pick first block on mesh
            auto tBlock = aMesh->get_block_by_index( 0 );

            // how many cells exist on current proc
            auto tNumberOfElements = tBlock->get_number_of_cells();

            // create nodes for these elements
            auto tNumberOfNodes = tBlock->get_number_of_vertices();

            // create node objects
            mNodes.resize(  tNumberOfNodes, nullptr );

             for( luint k=0; k<tNumberOfNodes; ++k )
             {
             mNodes( k ) = new fem::Node( tBlock->get_vertex_by_index( k ) );
             }

            // create equation objects
            mElements.resize( tNumberOfElements, nullptr );

            for( luint k=0; k<tNumberOfElements; ++k )
            {
                mElements( k ) = new fem::Element(
                        tBlock->get_cell_by_index( k ),
                        & aIWG,
                        mNodes,
                        aWeakBCs );

                mElements( k )->compute_jacobian_and_residual();
            }

            // create map for MSI
            map< moris_id, moris_index > tAdofMap;
            tBlock->get_adof_map( tAdofMap );

            // this part does not work yet in parallel
            auto tMSI = new moris::MSI::Model_Solver_Interface(
                    mElements,
                    aMesh->get_communication_table(),
                    tAdofMap,
                    tBlock->get_number_of_adofs_used_by_proc() );

            // create interface
            moris::MSI::MSI_Solver_Interface *  tSolverInput;
            tSolverInput = new moris::MSI::MSI_Solver_Interface( tMSI, tMSI->get_dof_manager() );

            // crete non-linear solver
            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Solver > tNonLinSolver = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            // crete linear solver
            moris::Solver_Factory  tSolFactory;

            // create solver object
            auto tLin = tSolFactory.create_solver( tSolverInput, SolverType::TRILINOSTEST );

            tNonLinSolver->set_linear_solver( tLin );

            // solve problem
            tNonLinSolver->solver_nonlinear_system();

            Matrix< DDRMat > tDOFs;

            tLin->import();
            tLin->get_solution_full( tDOFs );

            // write result into output
            //tLin->get_solution( tDOFs );

            uint tLength = tDOFs.length();

            // make sure that length of vector is correct
            MORIS_ASSERT( tLength == (uint) tBlock->get_number_of_adofs_used_by_proc(),
                    "Number of ADOFs does not match" );

            // fixme this is only temporary. Needed for integration error
            for( auto tElement : mElements )
            {
                 tElement->extract_values( tLin );
            }

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
