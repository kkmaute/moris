#include "cl_HMR_Lagrange_Node_Interpolation.hpp"

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

    void
    Lagrange_Node_Interpolation::set_weights( const Mat< real > & aTMatrix )
    {
        mWeights = aTMatrix;
    }

// ----------------------------------------------------------------------------

    const Mat< real > *
    Lagrange_Node_Interpolation::get_weights() const
    {
        return & mWeights;
    }

// ----------------------------------------------------------------------------

    void
    Lagrange_Node_Interpolation::set_coefficients(
            Cell< mtk::Vertex* > & aCoefficients )
    {
        mCoefficients = aCoefficients;
    }

// ----------------------------------------------------------------------------

    Cell< mtk::Vertex* > &
    Lagrange_Node_Interpolation::get_coefficients()
    {
        return mCoefficients;
    }

// ----------------------------------------------------------------------------

    const Cell< mtk::Vertex* > &
    Lagrange_Node_Interpolation::get_coefficients() const
    {
        return mCoefficients;
    }

// ----------------------------------------------------------------------------

    uint
    Lagrange_Node_Interpolation::get_number_of_coefficients() const
    {
        return mCoefficients.size();
    }

// ----------------------------------------------------------------------------

    Mat< moris_id >
    Lagrange_Node_Interpolation::get_ids() const
    {
        // get number of basis
        uint tNumberOfBasis = this->get_number_of_coefficients();

        // allocate output matrix
        Mat< moris_id > aIDs( tNumberOfBasis, 1 );

        // loop over all basis
        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            aIDs( k ) = mCoefficients( k )->get_id();
        }

        // return id matrix
        return aIDs;
    }

// ----------------------------------------------------------------------------

    Mat< moris_index >
    Lagrange_Node_Interpolation::get_indices() const
    {
        // get number of basis
        uint tNumberOfBasis = this->get_number_of_coefficients();

        // allocate output matrix
        Mat< moris_index > aIndices( tNumberOfBasis, 1 );

        // loop over all basis
        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            aIndices( k ) = mCoefficients( k )->get_index();
        }

        // return id matrix
        return aIndices;
    }

// ----------------------------------------------------------------------------

    Mat< uint >
    Lagrange_Node_Interpolation::get_owners() const
    {
        // get number of basis
        uint tNumberOfBasis = this->get_number_of_coefficients();

        // allocate output matrix
        Mat< uint > aOwners( tNumberOfBasis, 1 );

        // loop over all basis
        for( uint k=0; k<tNumberOfBasis; ++k )
        {
            aOwners( k ) = mCoefficients( k )->get_owner();
        }

        // return id matrix
        return aOwners;
    }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
