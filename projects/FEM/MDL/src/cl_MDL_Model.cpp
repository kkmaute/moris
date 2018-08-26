#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src
#include "cl_Solver_Factory.hpp"
#include "cl_Solver_Input.hpp"

#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Node_Obj.hpp"
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

            // create nodes for these elememnts
            auto tNumberOfNodes = tBlock->get_number_of_vertices();

            // FIXME : Check if finalize influences the indices
            Cell< MSI::Node* > tNodes( tNumberOfNodes, nullptr );
            for( luint k=0; k<tNumberOfNodes; ++k )
            {
                tNodes( k ) = new MSI::Node( tBlock->get_vertex_by_index( k ) );
            }

            // create equation objects
            Cell< MSI::Equation_Object* > tElements( tNumberOfElements, nullptr );

            for( luint k=0; k<tNumberOfElements; ++k )
            {
                tElements( k ) = new fem::Element(
                        tBlock->get_cell_by_index( k ),
                        & aIWG,
                        tNodes,
                        aWeakBCs );
            }

            // create interface object

            // this part does not work yet in parallel
            auto tMSI = new moris::MSI::Model_Solver_Interface(
                    tElements,
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

            // write result into output
            tLin->get_solution( aDOFs );

            // tidy up
            delete tSolverInput;

            // delete interface
            delete tMSI;

            // delete elements
            for( auto tElement : tElements )
            {
                delete tElement;
            }

            // delete nodes
            for( auto tNode : tNodes )
            {
                delete tNode;
            }

        }

//------------------------------------------------------------------------------

        Model::~Model()
        {


        }

//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */
