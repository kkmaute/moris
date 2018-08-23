#include "cl_MDL_Model.hpp"
#include "cl_FEM_Element.hpp"               //FEM/INT/src
namespace moris
{
    namespace mdl
    {
//------------------------------------------------------------------------------

        Model::Model(
                      mtk::Mesh   & aMesh,
                      fem::IWG    & aIWG,
                const Mat< real > & aInput,
                      Mat< real > & aResult ) :
                        mResult( aResult )
        {
            // pick first block on mesh
            auto tBlock = aMesh.get_block_by_index( 0 );

            // how many cells exist on current proc
            auto tNumberOfCells = tBlock->get_number_of_cells();

            // flag elements on this block
            for( luint e=0; e<tNumberOfCells; ++e )
            {

                // flag cell
                tBlock->get_cell_by_index( e )->set_t_matrix_flag();
            }

            // finalize mesh
            aMesh.finalize();

            // create nodes for these elememnts
            auto tNumberOfVertices = tBlock->get_number_of_vertices();

            // FIXME : Check if finalize influences the indices
            mNodes.resize( tNumberOfVertices, nullptr );

            for( luint k=0; k<tNumberOfVertices; ++k )
            {
                mNodes( k ) = new MSI::Node( tBlock->get_vertex_by_index( k ) );
            }

           // create elements from block
            mEquationObjects.resize( tNumberOfCells, nullptr );

            for( luint k=0; k<tNumberOfCells; ++k )
            {
                mEquationObjects( k ) = new fem::Element(
                        tBlock->get_cell_by_index( k ),
                        & aIWG,
                        mNodes,
                        aInput );
            }

            // solve model
            if( par_size() == 1)
            {
                // this part does not work yet in parallel
                moris::MSI::Model_Solver_Interface tMSI(
                        mEquationObjects,
                        aMesh.get_communication_table());
                tMSI.solve_system( mResult );
            }

            // write result in output
            // copy node values from equation object
            //for ( auto tElement : tListEqnObj )
            //{
            //    tElement->get_pdof_values( tNodeValues );
           // }
        }

//------------------------------------------------------------------------------

        Model::~Model()
        {

           // delete equation objects
           for( auto tObject : mEquationObjects )
           {
               delete tObject;
           }

           // delete nodes if they exist
           for( auto tNode : mNodes )
           {
               delete tNode;
           }
        }

//------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */
