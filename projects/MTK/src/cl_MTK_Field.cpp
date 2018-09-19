#include "fn_dot.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Block.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

    void
    Field::evaluate_node_values( const Matrix< DDRMat > & aCoefficients )
    {
        // ask mesh for number of nodes
        uint tNumberOfNodes = mBlock->get_number_of_vertices();

        Matrix< DDRMat > & tNodeValues = this->get_node_values();

        // allocate memory for matrix
        tNodeValues.set_size( tNumberOfNodes, mNumberOfDimensions );



        MORIS_ERROR( mNumberOfDimensions == 1,
                     "currently, only scalar fields are supported" );

        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            // get pointer to node
            auto tNode = mBlock->get_vertex_by_index( k );

            // get PDOFs from node
            auto tBSplines = tNode->get_interpolation()->get_coefficients();

            // get T-Matrix
            const Matrix< DDRMat > & tTMatrix = *tNode->get_interpolation()->get_weights();

            // get number of coefficients
            uint tNumberOfCoeffs = tTMatrix.length();

            // fill coeffs vector
            Matrix< DDRMat > tCoeffs( tNumberOfCoeffs, 1 );
            for( uint i=0; i<tNumberOfCoeffs; ++i )
            {
                tCoeffs( i ) = aCoefficients( tBSplines( i )->get_index() );
            }

            // write value into solution
            tNodeValues( k ) = dot( tTMatrix, tCoeffs );
        }
    }

//------------------------------------------------------------------------------

    void
    Field::evaluate_scalar_function(
                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
    {
        // get pointer to node values
        Matrix< DDRMat > & tNodeValues = this->get_node_values();

        // get number of nodes on block
        uint tNumberOfVertices = mBlock->get_number_of_vertices();

        // set size of node values
        tNodeValues.set_size( tNumberOfVertices, 1 );

        // loop over all vertices
        for( uint k=0; k<tNumberOfVertices; ++k )
        {
            // evaluate function at vertex cooridinates

            tNodeValues( k ) = aFunction(
                    mBlock->get_vertex_by_index( k )->get_coords() );

        }
    }

//------------------------------------------------------------------------------

    void
    Field::evaluate_node_values()
    {
        this->evaluate_node_values( mCoefficients );
    }

//------------------------------------------------------------------------------

    uint
    Field::get_interpolation_order() const
    {
        return this->mBlock->get_interpolation_order();
    }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
