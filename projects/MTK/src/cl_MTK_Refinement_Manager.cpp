#include "assert.hpp"
#include "cl_MTK_Refinement_Manager.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Block.hpp"
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        void
        Refinement_Manager::find_volume_and_surface_cells(
                Mat< moris_index > & aVolumeCellIndices,
                Mat< moris_index > & aSurfaceCellIndices,
          const Field              * aScalarField,
          const real                 aLowerBound,
          const real                 aUpperBound )
        {
            // this function only works if the field is a scalar field
            MORIS_ERROR( aScalarField->get_number_of_dimensions() == 1,
                    "find_intersected_elements() can only be done for a scalar field" );

            // get matrix to field values
            const Mat< real > & tVertexValues = * aScalarField->get_node_values();

            // get pointer to block on which field is defined
            auto tBlock = aScalarField->get_block();

            // make sure that node values are calculated
            MORIS_ERROR( tVertexValues.length() == tBlock->get_number_of_vertices(),
                    "number of field values does not match number of vertices on block" );

            // get number of cells on this block
            uint tNumberOfCells = tBlock->get_number_of_cells();

            std::cout << "Cells: " << tNumberOfCells << std::endl;

            // get number of vertices per cell ( must be the same for all elements on block )
            uint tNumberOfVerticesPerCell
                = tBlock->get_cell_by_index( 0 )->get_number_of_vertices();

            // allocate matrix with values
            Mat< real > tCellValues( tNumberOfVerticesPerCell, 1 );

            // counter for cells
            uint tSurfaceCount = 0;
            uint tVolumeCount = 0;

            // allocate vector for output elements
            aVolumeCellIndices.set_size( tNumberOfCells, 1 );
            aSurfaceCellIndices.set_size( tNumberOfCells, 1 );

            // loop over all cells on this block
            for( uint e=0; e<tNumberOfCells; ++e )
            {
                // get pointer to cell
                const Cell * tCell = tBlock->get_cell_by_index( e );

                // get cell of vertex pointers
                Mat< moris_index > tVertices = tCell->get_vertex_indices();

                // copy values of associated vertices into field
                for( uint k=0; k<tNumberOfVerticesPerCell; ++k )
                {
                    tCellValues( k ) = tVertexValues( tVertices( k ) );
                }

                // test if cell is internal
                if( tCellValues.max() <= aUpperBound )
                {
                    aVolumeCellIndices( tVolumeCount++ ) = e;

                    // jump to next cell
                    continue;
                }

                // test if cell is intersected
                if ( tCellValues.min() <= aLowerBound && tCellValues.max() >= aUpperBound )
                {
                    aSurfaceCellIndices( tSurfaceCount++ ) = e;
                }
            }

            // resize output vector
            aVolumeCellIndices.resize( tVolumeCount, 1 );
            aSurfaceCellIndices.resize( tSurfaceCount, 1 );
        }

//------------------------------------------------------------------------------
    }
}
