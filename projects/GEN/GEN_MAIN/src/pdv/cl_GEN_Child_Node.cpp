#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell_Info.hpp"

namespace moris
{
namespace ge
{

    //--------------------------------------------------------------------------------------------------------------

    Child_Node::Child_Node(
        Matrix< DDUMat >           aAncestorNodeIndices,
        Cell< Matrix< DDRMat > >   aAncestorNodeCoordinates,
        const xtk::Basis_Function& aBasisFunction,
        Matrix< DDRMat >           aLocalCoordinatesInAncestor ) :
        mAncestorNodeIndices( aAncestorNodeIndices ),
        mAncestorNodeCoordinates( aAncestorNodeCoordinates )
    {
        if ( aAncestorNodeIndices.length() == 0 )
        {
            MORIS_ERROR( 0, "BUG" );
        }
        // Check that ancestor info is consistent
        MORIS_ASSERT( aAncestorNodeIndices.length() == aAncestorNodeCoordinates.size(),
            "Number of ancestor node indices must be consistent with the number of sets of node coordinates" );

        // Evaluate basis function
        aBasisFunction.evaluate_basis_function( aLocalCoordinatesInAncestor, mBasisValues );
    }

    Child_Node::Child_Node( moris::mtk::Cell* aCell,
        Matrix< DDRMat >*                     aLocalCoordinates )
    {
        mParentCell       = aCell;
        mLocalCoordinates = *aLocalCoordinates;

        Cell< moris::mtk::Vertex* > tVerts = mParentCell->get_vertex_pointers();

        mAncestorNodeIndices.resize( 1, tVerts.size() );
        mAncestorNodeCoordinates.resize( tVerts.size(), Matrix< DDRMat >( 1, 3 ) );

        for ( moris::uint iCast = 0; iCast < tVerts.size(); iCast++ )
        {
            mAncestorNodeIndices( iCast )     = (moris::uint)tVerts( iCast )->get_index();
            mAncestorNodeCoordinates( iCast ) = tVerts( iCast )->get_coords();
        }

        aCell->get_cell_info()->eval_N( mLocalCoordinates, mBasisValues );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Child_Node::interpolate_field_value( Field* aField )
    {
        // Get field values from ancestors
        Matrix< DDRMat > tGeometryFieldValues( mAncestorNodeIndices.length(), 1 );
        for ( uint tAncestorNode = 0; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++ )
        {
            tGeometryFieldValues( tAncestorNode ) = aField->get_field_value( mAncestorNodeIndices( tAncestorNode ), mAncestorNodeCoordinates( tAncestorNode ) );
        }

        // Return interpolated value
        return Matrix< DDRMat >( mBasisValues * tGeometryFieldValues )( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Child_Node::join_field_sensitivities( Field* aField )
    {
        // Initialize using first ancestor
        mJoinedSensitivities = aField->get_dfield_dadvs(
            mAncestorNodeIndices( 0 ),
            mAncestorNodeCoordinates( 0 ) );
        mJoinedSensitivities = mJoinedSensitivities * mBasisValues( 0 );

        // Get sensitivity values from other ancestors
        for ( uint tAncestorNode = 1; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++ )
        {
            // Get scaled sensitivities
            const Matrix< DDRMat >& tAncestorSensitivities = mBasisValues( tAncestorNode ) * aField->get_dfield_dadvs( mAncestorNodeIndices( tAncestorNode ), mAncestorNodeCoordinates( tAncestorNode ) );

            // Join sensitivities
            uint tJoinedSensitivityLength = mJoinedSensitivities.n_cols();
            mJoinedSensitivities.resize( 1, tJoinedSensitivityLength + tAncestorSensitivities.length() );
            for ( uint tAncestorSensitivity = 0; tAncestorSensitivity < tAncestorSensitivities.length(); tAncestorSensitivity++ )
            {
                mJoinedSensitivities( tJoinedSensitivityLength + tAncestorSensitivity ) = tAncestorSensitivities( tAncestorSensitivity );
            }
        }

        return mJoinedSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    Child_Node::join_determining_adv_ids( Field* aField )
    {
        // Initialize using first ancestor
        Matrix< DDSMat > tJoinedDeterminingADVs = aField->get_determining_adv_ids(
            mAncestorNodeIndices( 0 ),
            mAncestorNodeCoordinates( 0 ) );

        // Get sensitivity values from other ancestors
        for ( uint tAncestorNode = 1; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++ )
        {
            // Get scaled sensitivities
            Matrix< DDSMat > tAncestorDependingADVs = aField->get_determining_adv_ids(
                mAncestorNodeIndices( tAncestorNode ),
                mAncestorNodeCoordinates( tAncestorNode ) );
            // Join sensitivities
            uint tJoinedADVLength = tJoinedDeterminingADVs.n_cols();
            tJoinedDeterminingADVs.resize( 1, tJoinedADVLength + tAncestorDependingADVs.length() );
            for ( uint tAncestorADV = 0; tAncestorADV < tAncestorDependingADVs.length(); tAncestorADV++ )
            {
                tJoinedDeterminingADVs( tJoinedADVLength + tAncestorADV ) = tAncestorDependingADVs( tAncestorADV );
            }
        }

        return tJoinedDeterminingADVs;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Child_Node::get_dfield_dcoordinates(
        Field*            aField,
        Matrix< DDRMat >& aSensitivities )
    {
        MORIS_ERROR( false, "A child node that isn't an intersection node does not know dfield_dcoordinates" );
    }

    //--------------------------------------------------------------------------------------------------------------

}// namespace ge
}// namespace moris
