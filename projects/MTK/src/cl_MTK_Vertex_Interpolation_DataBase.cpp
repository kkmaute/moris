/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Interpolation_DataBase.cpp
 *
 */

#include "cl_MTK_Vertex_Interpolation_DataBase.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "fn_TOL_Capacities.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        Vertex_Interpolation_DataBase::Vertex_Interpolation_DataBase( moris_index aVertexIndex,
            moris_index                                                           aVertexOrder,
            mtk::Mesh*                                                            aMesh ) :
            mVertexIndex( aVertexIndex ),
            mVertexOrder( aVertexOrder ),
            mMesh( aMesh )
        {
            this->set_outward_data();
        }

        //------------------------------------------------------------------------------

        Matrix< IdMat >
        Vertex_Interpolation_DataBase::get_ids() const
        {
            // get the basis pointer and number of elements
            moris_id*   tBasisIdPtr = mMesh->get_basis_ids( mVertexIndex, mVertexOrder );
            moris_index tOffset     = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

            Matrix< IdMat > tBasisIds( tBasisIdPtr, tOffset, 1, false, true );

            return tBasisIds;
        }

        //------------------------------------------------------------------------------

        Matrix< IndexMat >
        Vertex_Interpolation_DataBase::get_indices() const
        {
            // get the basis pointer and number of elements
            moris_index* tBasisIndicesPtr = mMesh->get_basis_indicies( mVertexIndex, mVertexOrder );
            moris_index  tOffset          = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

            Matrix< IdMat > tBasisIndices( tBasisIndicesPtr, tOffset, 1, false, true );

            return tBasisIndices;
        }

        //------------------------------------------------------------------------------

        Matrix< IdMat >
        Vertex_Interpolation_DataBase::get_owners() const
        {
            // get the basis pointer and number of elements
            moris_id*   tBasisOwnersPtr = mMesh->get_basis_owners( mVertexIndex, mVertexOrder );
            moris_index tOffset         = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

            Matrix< IdMat > tBasisOwners( tBasisOwnersPtr, tOffset, 1, false, true );

            return tBasisOwners;
        }

        //------------------------------------------------------------------------------

        void
        Vertex_Interpolation_DataBase::set_weights( const moris::Matrix< DDRMat >& aWeights )
        {
            MORIS_ERROR( 0, "set_weights not implemented in fem vertex interpolation" );
            return;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat >*
        Vertex_Interpolation_DataBase::get_weights() const
        {
            return &mWeights;
        }

        //------------------------------------------------------------------------------

        void
        Vertex_Interpolation_DataBase::set_coefficients( Vector< Vertex* >& aCoefficients )
        {
            MORIS_ERROR( 0, "set_coefficients not implemented in fem vertex interpolation" );
            return;
        }

        //------------------------------------------------------------------------------

        Vector< Vertex* >&
        Vertex_Interpolation_DataBase::get_coefficients()
        {
            MORIS_ERROR( 0, "get_coefficients not implemented in fem vertex interpolation" );
            return mCoefficients;
        }

        //------------------------------------------------------------------------------

        const Vector< Vertex* >&
        Vertex_Interpolation_DataBase::get_coefficients() const
        {
            MORIS_ERROR( 0, "get_coefficients not implemented in fem vertex interpolation" );
            return mCoefficients;
        }

        //------------------------------------------------------------------------------

        uint
        Vertex_Interpolation_DataBase::get_number_of_coefficients() const
        {
            MORIS_ERROR( 0, "get_number_of_coefficients not implemented in fem vertex interpolation" );
            return 0;
        }

        //------------------------------------------------------------------------------

        void
        Vertex_Interpolation_DataBase::set_outward_data()
        {
            moris::real* tBasisWightsPtr = mMesh->get_basis_weights( mVertexIndex, mVertexOrder );
            moris_index  tOffset         = mMesh->get_basis_length( mVertexIndex, mVertexOrder );

            mWeights = Matrix< DDRMat >( tBasisWightsPtr, tOffset, 1, false, true );
        }

        //------------------------------------------------------------------------------

        size_t
        Vertex_Interpolation_DataBase::capacity()
        {
            size_t tCapacity = 0;

            tCapacity += sizeof( mVertexIndex );
            tCapacity += sizeof( mVertexOrder );
            tCapacity += sizeof( mMesh );

            tCapacity += sizeof( mWeights ) + mWeights.capacity();

            return tCapacity;
        }
    }// namespace mtk
}// namespace moris
