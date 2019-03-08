#include "../cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation_STK.hpp"

namespace moris
{
    namespace mtk
    {
// ----------------------------------------------------------------------------

    void
    Node_Interpolation_STK::set_weights( const Matrix< DDRMat > & aTMatrix )
    {
        mWeights = aTMatrix;
    }

// ----------------------------------------------------------------------------

    const Matrix< DDRMat > *
    Node_Interpolation_STK::get_weights() const
    {
        return & mWeights;
    }

// ----------------------------------------------------------------------------

    void
    Node_Interpolation_STK::set_coefficients( moris::Cell< mtk::Vertex* > & aCoefficients )
    {
        mCoefficients = aCoefficients;
    }

// ----------------------------------------------------------------------------

    moris::Cell< mtk::Vertex* > &
    Node_Interpolation_STK::get_coefficients()
    {
        return mCoefficients;
    }

// ----------------------------------------------------------------------------

    const moris::Cell< mtk::Vertex* > &
    Node_Interpolation_STK::get_coefficients() const
    {
        return mCoefficients;
    }

// ----------------------------------------------------------------------------

    uint
    Node_Interpolation_STK::get_number_of_coefficients() const
    {
        return mCoefficients.size();
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
            aIDs( k ) =  mVertex->get_id();
        }

        // return id matrix
        return aIDs;
    }

// ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Node_Interpolation_STK::get_indices() const
    {
        // get number of basis
        uint tNumberOfBasis = this->get_number_of_coefficients();

        // allocate output matrix
        Matrix< IndexMat > aIndices( tNumberOfBasis, 1 );

        // loop over all basis
        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            aIndices( k ) = mVertex->get_index();
        }

        // return id matrix
        return aIndices;
    }

    Matrix< IdMat >
    Node_Interpolation_STK::get_owners() const
    {
        // get number of basis
        uint tNumberOfBasis = this->get_number_of_coefficients();

        // allocate output matrix
        Matrix< IdMat > aOwners( tNumberOfBasis, 1 );

        // loop over all basis
        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            aOwners( k ) = mVertex->get_owner();
        }

        // return id matrix
        return aOwners;
    }
    } /* namespace hmr */
} /* namespace moris */
