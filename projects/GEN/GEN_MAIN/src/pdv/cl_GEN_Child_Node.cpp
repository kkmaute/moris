#include "cl_GEN_Child_Node.hpp"
#include "cl_GEN_Field.hpp"

#include "cl_MTK_Cell.hpp"                              //MTK/src
#include "cl_MTK_Vertex.hpp"                            //MTK/src
#include "cl_MTK_Interpolation_Function_Base.hpp"       //MTK/src
#include "cl_MTK_Interpolation_Function_Factory.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"                             //MTK/src
#include "cl_MTK_Cell_Info.hpp"                         //MTK/src
#include "cl_MTK_Cell_Info_Quad4.hpp"                   //MTK/src
#include "cl_MTK_Cell_Info_Hex8.hpp"                    //MTK/src

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Child_Node::Child_Node(
                Matrix< DDUMat >          aAncestorNodeIndices,
                Cell< Matrix< DDRMat > >  aAncestorNodeCoordinates,
                Element_Intersection_Type aBasisFunction,
                Matrix< DDRMat >          aLocalCoordinatesInAncestor )
                : mAncestorNodeIndices( aAncestorNodeIndices )
                , mAncestorNodeCoordinates( aAncestorNodeCoordinates )
        {
            if ( aAncestorNodeIndices.length() == 0 )
            {
                MORIS_ERROR( 0, "BUG" );
            }
            // Check that ancestor info is consistent
            MORIS_ASSERT( aAncestorNodeIndices.length() == aAncestorNodeCoordinates.size(),
                    "Number of ancestor node indices must be consistent with the number of sets of node coordinates" );

            // construct interpolator and evaluate basis function
            mtk::Interpolation_Function_Factory tFactory;

            mtk::Interpolation_Function_Base* tInterpolation;

            switch ( aBasisFunction )
            {
                case Element_Intersection_Type::Linear_1D:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::LINE,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tInterpolation->eval_N( aLocalCoordinatesInAncestor, mBasisValues );
                    break;
                }
                case Element_Intersection_Type::Linear_2D:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tInterpolation->eval_N( aLocalCoordinatesInAncestor, mBasisValues );
                    break;
                }
                case Element_Intersection_Type::Linear_3D:
                {
                    tInterpolation = tFactory.create_interpolation_function(
                            mtk::Geometry_Type::HEX,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    tInterpolation->eval_N( aLocalCoordinatesInAncestor, mBasisValues );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Child_Node::Child_Node - Interpolation type not implemented." );
                }
            }

            // delete interpolator
            delete tInterpolation;
        }

        //--------------------------------------------------------------------------------------------------------------

        Child_Node::Child_Node( moris::mtk::Cell* aCell,
                Matrix< DDRMat >*                 aLocalCoordinates,
                bool                              aEvaluateAsLinear )
        {
            mParentCell       = aCell;
            mLocalCoordinates = *aLocalCoordinates;

            // use specifically bi or tri-linear interpolation to compute level set value
            if ( aEvaluateAsLinear )
            {
                Cell< moris::mtk::Vertex* > tVerts = mParentCell->get_vertex_pointers();

                if ( mParentCell->get_cell_info()->get_cell_geometry() == mtk::Geometry_Type::QUAD )
                {
                    mAncestorNodeIndices.resize( 1, 4 );
                    mAncestorNodeCoordinates.resize( 4, Matrix< DDRMat >( 1, 2 ) );

                    for ( moris::uint iCast = 0; iCast < 4; iCast++ )
                    {
                        mAncestorNodeIndices( iCast )     = (moris::uint)tVerts( iCast )->get_index();
                        mAncestorNodeCoordinates( iCast ) = tVerts( iCast )->get_coords();
                    }

                    mtk::Cell_Info_Quad4 tQuad4CellInfo;
                    tQuad4CellInfo.eval_N( mLocalCoordinates, mBasisValues );
                }

                if ( mParentCell->get_cell_info()->get_cell_geometry() == mtk::Geometry_Type::HEX )
                {
                    mAncestorNodeIndices.resize( 1, 8 );
                    mAncestorNodeCoordinates.resize( 8, Matrix< DDRMat >( 1, 3 ) );

                    for ( moris::uint iCast = 0; iCast < 8; iCast++ )
                    {
                        mAncestorNodeIndices( iCast )     = (moris::uint)tVerts( iCast )->get_index();
                        mAncestorNodeCoordinates( iCast ) = tVerts( iCast )->get_coords();
                    }

                    mtk::Cell_Info_Hex8 tLinearCellInfo;
                    tLinearCellInfo.eval_N( mLocalCoordinates, mBasisValues );
                }
            }
            else
            {
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
        }

        //--------------------------------------------------------------------------------------------------------------

        real
        Child_Node::interpolate_field_value( Field* aField )
        {
            // Get field values from ancestors
            Matrix< DDRMat > tGeometryFieldValues( mAncestorNodeIndices.length(), 1 );
            for ( uint tAncestorNode = 0; tAncestorNode < mAncestorNodeIndices.length(); tAncestorNode++ )
            {
                tGeometryFieldValues( tAncestorNode ) =
                        aField->get_field_value(
                                mAncestorNodeIndices( tAncestorNode ),
                                mAncestorNodeCoordinates( tAncestorNode ) );
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
                const Matrix< DDRMat > tAncestorSensitivities =
                        mBasisValues( tAncestorNode ) * aField->get_dfield_dadvs( mAncestorNodeIndices( tAncestorNode ), mAncestorNodeCoordinates( tAncestorNode ) );

                // Join sensitivities
                uint tJoinedSensitivityLength = mJoinedSensitivities.n_cols();

                // FIXME: excessive resize
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

    }    // namespace ge
}    // namespace moris
