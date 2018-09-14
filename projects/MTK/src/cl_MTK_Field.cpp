#include "fn_dot.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Block.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

    void
    Field::evaluate_node_values( const Mat< real > & aCoefficients )
    {
        // ask mesh for number of nodes
        uint tNumberOfNodes = mBlock->get_number_of_vertices();

        // allocate memory for matrix
        mNodeValues->set_size( tNumberOfNodes, mNumberOfDimensions );

        Mat< real > & tNodeValues = * mNodeValues;

        MORIS_ERROR( mNumberOfDimensions == 1,
                     "currently, only scalar fields are supported" );

        for( uint k=0; k<tNumberOfNodes; ++k )
        {
            // get pointer to node
            auto tNode = mBlock->get_vertex_by_index( k );

            // get PDOFs from node
            auto tBSplines = tNode->get_adof_pointers();

            // get T-Matrix
            const Mat< real > & tTMatrix = *tNode->get_t_matrix();

            // get number of coefficients
            uint tNumberOfCoeffs = tTMatrix.length();

            // fill coeffs vector
            Mat< real > tCoeffs( tNumberOfCoeffs, 1 );
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
    Field::evaluate_node_values()
    {
        this->evaluate_node_values( * mCoefficients );
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
