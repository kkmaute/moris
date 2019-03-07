#include "cl_MTK_Node_Interpolation_STK.hpp"
#include "../cl_MTK_Vertex.hpp"

namespace moris
{
    namespace mtk
    {
// ----------------------------------------------------------------------------

    void
    Node_Interpolation_STK::set_weights( const Matrix< DDRMat > & aTMatrix )
    {
        mWeightsParent = aTMatrix;
    }

// ----------------------------------------------------------------------------

    const Matrix< DDRMat > *
	Node_Interpolation_STK::get_weights() const
    {
        return & mWeightsParent;
    }

// ----------------------------------------------------------------------------

    void
    Node_Interpolation_STK::set_coefficients( Cell< mtk::Vertex* > & aCoefficients )
    {
        mCoefficientsParent = aCoefficients;
    }

// ----------------------------------------------------------------------------

    Cell< mtk::Vertex* > &
    Node_Interpolation_STK::get_coefficients()
    {
        return mCoefficientsParent;
    }

// ----------------------------------------------------------------------------

    const Cell< mtk::Vertex* > &
    Node_Interpolation_STK::get_coefficients() const
    {
        return mCoefficientsParent;
    }

// ----------------------------------------------------------------------------

    uint
    Node_Interpolation_STK::get_number_of_coefficients() const
    {
        return mCoefficientsParent.size();
    }

// ----------------------------------------------------------------------------

    Matrix< IdMat >
    Node_Interpolation_STK::get_ids() const
    {
        // get number of basis
        uint tNumberOfBasis = this->get_number_of_coefficients();

        // allocate output matrix
        Matrix< IdMat > aIDs( tNumberOfBasis, 1 );

        // loop over all basis
        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            aIDs( k ) = mCoefficientsParent( k )->get_id();
        }

        // return id matrix
        return aIDs;
    }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
