#include "cl_Stopwatch.hpp" //CHR/src

#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src
#include "cl_FEM_Node.hpp"               //FEM/INT/src

#include "cl_Solver_Factory.hpp"
#include "cl_Solver_Input.hpp"

#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
//#include "cl_MSI_Node_Obj.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"

namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        Model::Model(
                mtk::Mesh         & aMesh,
                fem::IWG          & aIWG,
                const Mat< real > & aWeakBCs,
                Mat< real >       & aDOFs )
        {


            // pick first block on mesh
            auto tBlock = aMesh.get_block_by_index( 0 );

            // how many cells exist on current proc
            auto tNumberOfElements = tBlock->get_number_of_cells();

            // flag elements on this block
            for( luint e=0; e<tNumberOfElements; ++e )
            {
                // flag cell
                tBlock->get_cell_by_index( e )->set_t_matrix_flag();
            }

            // finalize mesh ( calculate T-Matrices etc )
            aMesh.finalize();

            // create nodes for these elements
            auto tNumberOfNodes = tBlock->get_number_of_vertices();

            std::cout << " Model " << tNumberOfElements << " " << tNumberOfNodes << std::endl;

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

            // create interface object

            // this part does not work yet in parallel
            auto tMSI = new moris::MSI::Model_Solver_Interface(
                    mElements,
                    aMesh.get_communication_table() );

            // create interface
            moris::MSI::MSI_Solver_Interface *  tSolverInput;
            tSolverInput = new moris::MSI::MSI_Solver_Interface( tMSI, tMSI->get_dof_manager() );

            // crete linear solver
            moris::Solver_Factory  tSolFactory;

            // create solver object
            auto tLin = tSolFactory.create_solver( tSolverInput );


            // solve problem
            tLin->solve_linear_system();


            // fixme this is only temporary. Needed for integration error
            for( auto tElement : mElements )
            {
                tElement->get_pdof_values( tLin );
            }


            // write result into output
            tLin->get_solution( aDOFs );

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
                real (*aFunction)( const Mat< real > & aPoint ) )
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
