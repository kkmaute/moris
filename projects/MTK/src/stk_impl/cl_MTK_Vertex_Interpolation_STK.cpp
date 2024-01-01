/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_STK.cpp
 *
 */

#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation_STK.hpp"

namespace moris
{
    namespace mtk
    {
// ----------------------------------------------------------------------------

    void
	Vertex_Interpolation_STK::set_weights( const Matrix< DDRMat > & aTMatrix )
    {
        mWeights = aTMatrix;
    }

// ----------------------------------------------------------------------------

    const Matrix< DDRMat > *
	Vertex_Interpolation_STK::get_weights() const
    {
        return & mWeights;
    }

// ----------------------------------------------------------------------------

    void
	Vertex_Interpolation_STK::set_coefficients( moris::Vector< mtk::Vertex* > & aCoefficients )
    {
        mCoefficients = aCoefficients;
    }

// ----------------------------------------------------------------------------

    moris::Vector< mtk::Vertex* > &
	Vertex_Interpolation_STK::get_coefficients()
    {
        return mCoefficients;
    }

// ----------------------------------------------------------------------------

    const moris::Vector< mtk::Vertex* > &
	Vertex_Interpolation_STK::get_coefficients() const
    {
        return mCoefficients;
    }

// ----------------------------------------------------------------------------

    uint
	Vertex_Interpolation_STK::get_number_of_coefficients() const
    {
        //return mCoefficients.size();
        return 1;
    }

// ----------------------------------------------------------------------------

    Matrix< IdMat >
    Vertex_Interpolation_STK::get_ids() const
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
    Vertex_Interpolation_STK::get_indices() const
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
    Vertex_Interpolation_STK::get_owners() const
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

