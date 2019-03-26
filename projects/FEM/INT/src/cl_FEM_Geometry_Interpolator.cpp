#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "op_div.hpp"


#include "cl_FEM_Geometry_Interpolator.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Geometry_Interpolator::Geometry_Interpolator( const Interpolation_Rule & aInterpolationRule )
        {
            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // create member pointer to time interpolation function
            mTimeInterpolation  = aInterpolationRule.create_time_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceBases = mSpaceInterpolation->get_number_of_bases();
            mNumSpaceDim   = mSpaceInterpolation->get_number_of_dimensions();

            // number of time bases and dimensions
            mNumTimeBases = mTimeInterpolation->get_number_of_bases();
            mNumTimeDim   = mTimeInterpolation->get_number_of_dimensions();

            // set default xHat and tHat
            mXHat.set_size( mNumSpaceBases, mNumSpaceDim, 0.0);
            mTHat.set_size( mNumTimeBases,  mNumTimeDim,  0.0);

            // set member geometry type
            mGeometryType = aInterpolationRule.get_geometry_type();

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

        Geometry_Interpolator::Geometry_Interpolator( const Interpolation_Rule & aInterpolationRule,
                                                      const bool aSpaceSideset )
        {
        	// set bool for side interpolation to true
        	mSpaceSideset = true;

            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // create member pointer to time interpolation function
            mTimeInterpolation  = aInterpolationRule.create_time_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceBases = mSpaceInterpolation->get_number_of_bases();
            mNumSpaceDim   = mSpaceInterpolation->get_number_of_dimensions();

            // number of time bases and dimensions
            mNumTimeBases = mTimeInterpolation->get_number_of_bases();
            mNumTimeDim   = mTimeInterpolation->get_number_of_dimensions();

            // set default xHat and tHat
            mXHat.set_size( mNumSpaceBases, mNumSpaceDim, 0.0 );
            mTHat.set_size( mNumTimeBases,  mNumTimeDim,  0.0 );

            // set member geometry type
            mGeometryType = aInterpolationRule.get_geometry_type();

            // set member side geometry type
            this->get_auto_side_geometry_type();

            // create side interpolation rule
            Interpolation_Rule tSideInterpolationRule( mSideGeometryType,
                                                       aInterpolationRule.get_space_interpolation_type(),
                                                       aInterpolationRule.get_space_interpolation_order(),
                                                       aInterpolationRule.get_time_interpolation_type(),
                                                       aInterpolationRule.get_time_interpolation_order() );

            // create member pointer to side space interpolation function
            mSideSpaceInterpolation = tSideInterpolationRule.create_space_interpolation_function();

            // set pointers for second derivative depending on space and time dimensions
            this->set_function_pointers();
        }

//------------------------------------------------------------------------------

        Geometry_Interpolator::~Geometry_Interpolator()
        {
            // delete interpolation functions
            if( mSpaceInterpolation != NULL )
             {
                 delete mSpaceInterpolation;
             }

             if( mTimeInterpolation != NULL )
             {
                 delete mTimeInterpolation;
             }

             if( mSideSpaceInterpolation != NULL )
             {
                 delete mSideSpaceInterpolation;
             }
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::set_coeff( const Matrix< DDRMat > & aXHat,
                                               const Matrix< DDRMat > & aTHat )
        {
            //check the space coefficients input size
            MORIS_ASSERT( ( ( aXHat.n_cols() == mNumSpaceDim ) && ( aXHat.n_rows() == mNumSpaceBases )),
                          " Geometry_Interpolator::set_coeff - Wrong input size (aXHat). ");

            // set the space coefficients
            mXHat = aXHat;

            //check the time coefficients input size
            MORIS_ASSERT( ( ( aTHat.n_cols() == mNumTimeDim ) && ( aTHat.n_rows() == mNumTimeBases )),
                           " Geometry_Interpolator::set_coeff - Wrong input size (aTHat). ");

            // set the time coefficients
            mTHat = aTHat;
        }


//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::get_space_sideset_param_coords( const moris_index aSpaceOrdinal )
        {
            //FIXME check spaceOrdinal

            // get space sideset vertex indices
            moris::Cell< moris::moris_index > tVerticesIndices;
            switch ( mGeometryType )
            {
                case ( mtk::Geometry_Type::LINE ): //FIXME LINE2 only
                {
                    moris::Cell< moris::Cell< moris::moris_index > > tVerticesIndicesLine2
                        = { { 0 },
                            { 1 } };
                    tVerticesIndices.resize( 1, -1 );
                    tVerticesIndices = tVerticesIndicesLine2( aSpaceOrdinal );
                    break;
                }

                case ( mtk::Geometry_Type::QUAD ): //FIXME QUAD4 only
                {
                    moris::Cell< moris::Cell< moris::moris_index > > tVerticesIndicesQuad4
                        = { { 0, 1 },
                            { 1, 2 },
                            { 2, 3 },
                            { 3, 0 } };
                    tVerticesIndices.resize( 2, -1 );
                    tVerticesIndices = tVerticesIndicesQuad4( aSpaceOrdinal );
                    break;
                }

                case ( mtk::Geometry_Type::HEX ): // FIXME HEX8 only
                {
                    moris::Cell< moris::Cell< moris::moris_index > > tVerticesIndicesHex8
                        = { { 0, 1, 5, 4 },
                            { 1, 2, 6, 5 },
                            { 2, 3, 7, 6 },
                            { 0, 4, 7, 3 },
                            { 0, 3, 2, 1 },
                            { 4, 5, 6, 7 } };
                    tVerticesIndices.resize( 4, -1 );
                    tVerticesIndices = tVerticesIndicesHex8( aSpaceOrdinal );
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "Geometry_Interpolator::get_space_sideset_param_coords - undefined geometry type " );
                    break;
                }
            }

            // initialize the time sideset parametric coordinates matrix
            uint tNumOfVertices = tVerticesIndices.size();
            Matrix< DDRMat > tSpaceParamCoords = mSpaceInterpolation->get_param_coords();
            Matrix< DDRMat > tSidesetSpaceParamCoords( mNumSpaceDim, tNumOfVertices );
            for( uint i = 0; i < tNumOfVertices; i++ )
            {
                moris_index tTreatedVertex = tVerticesIndices( i );
                tSidesetSpaceParamCoords( { 0, mNumSpaceDim - 1 }, { i, i } )
                    = tSpaceParamCoords( { 0, mNumSpaceDim - 1 }, { tTreatedVertex, tTreatedVertex } );
            }

            // get a vector of ones
            Matrix< DDRMat > tParamCoords( mNumSpaceDim + mNumTimeDim, tNumOfVertices * mNumTimeBases );
            Matrix< DDRMat > tTimeParamCoords  = mTimeInterpolation->get_param_coords();
            Matrix< DDRMat > tOnes( 1, tNumOfVertices, 1.0 );
            for( uint i = 0; i < mNumTimeBases; i++ )
            {
                tParamCoords( {0, mNumSpaceDim - 1 }, { i*tNumOfVertices, (i+1)*tNumOfVertices-1 } )
                    = tSidesetSpaceParamCoords.matrix_data();

                // fill the space time parametric coordinates matrix with time coordinates
                tParamCoords( { mNumSpaceDim, mNumSpaceDim }, { i*tNumOfVertices, (i+1)*tNumOfVertices-1 } )
                    = tTimeParamCoords( i ) * tOnes;
            }

            // return the parametric coordinates of the space sideset
            return tParamCoords;
        }

//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::get_time_sideset_param_coords( const moris_index aTimeOrdinal )
        {
            //FIXME check on the time ordinal only, 0 or 1

            // initialize the time sideset parametric coordinates matrix
            Matrix< DDRMat > tParamCoords( mNumSpaceDim + mNumTimeDim, mNumSpaceBases );

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, mNumSpaceBases, 1.0 );

            // fill the space time parametric coordinates matrix with space coordinates
            tParamCoords( { 0, mNumSpaceDim-1 }, { 0, mNumSpaceBases-1 } )
               = mSpaceInterpolation->get_param_coords().matrix_data();

            // fill the space time parametric coordinates matrix with time coordinates
            tParamCoords( { mNumSpaceDim, mNumSpaceDim }, { 0, mNumSpaceBases-1 })
                = mTimeInterpolation->get_param_coords()( aTimeOrdinal ) * tOnes;

            // return parametric coordinates
            return tParamCoords;
        }

//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::get_space_time_param_coords()
        {
            // initialize the space time parametric coordinates matrix
            Matrix< DDRMat > tParamCoords( mNumSpaceDim + mNumTimeDim, mNumSpaceBases * mNumTimeBases );

            // get the soace parametric coordinates
            Matrix< DDRMat > tSpaceParamCoords = mSpaceInterpolation->get_param_coords();

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords  = mTimeInterpolation->get_param_coords();

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, mNumSpaceBases, 1.0 );

            // loop on the time bases
            for( uint i = 0; i < mNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tParamCoords( { 0, mNumSpaceDim-1 }, { i * mNumSpaceBases, ( i + 1 ) * mNumSpaceBases-1 })
                    = tSpaceParamCoords.matrix_data();

                // fill the space time parametric coordinates matrix with time coordinates
                tParamCoords( { mNumSpaceDim, mNumSpaceDim }, { i * mNumSpaceBases, ( i + 1 ) * mNumSpaceBases-1 })
                    = tTimeParamCoords( i )*tOnes;
            }

            // return the space time parametric coordinates
            return tParamCoords;
        }

//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::NXi( const Matrix< DDRMat > & aXi ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> tN = mSpaceInterpolation->eval_N( aXi );
            return tN;
         }

//------------------------------------------------------------------------------

         Matrix < DDRMat > Geometry_Interpolator::NTau( const Matrix< DDRMat > & aTau ) const
         {
             // pass data through interpolation function
             Matrix <DDRMat> tN = mTimeInterpolation->eval_N( aTau );
             return tN;
         }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::dNdXi( const Matrix< DDRMat > & aXi ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> tdNdXi = mSpaceInterpolation->eval_dNdXi( aXi );
            return tdNdXi;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::dNdTau( const Matrix< DDRMat > & aTau) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> tdNdTau = mTimeInterpolation->eval_dNdXi( aTau );
            return tdNdTau;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::d2NdXi2( const Matrix< DDRMat > & aXi ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> td2NdXi2 = mSpaceInterpolation->eval_d2NdXi2( aXi );
            return td2NdXi2;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::d2NdTau2( const Matrix< DDRMat > & aTau ) const
        {
            // pass data through interpolation function
            Matrix <DDRMat> td2NdTau2 = mTimeInterpolation->eval_d2NdXi2( aTau );
            return td2NdTau2;
        }
//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::space_jacobian( const Matrix< DDRMat > & adNdXi ) const
        {
            Matrix< DDRMat > tJt = adNdXi * mXHat ;
            return tJt;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::time_jacobian( const Matrix< DDRMat > & adNdTau ) const
        {
            Matrix< DDRMat > tJt = adNdTau * mTHat ;
            return tJt;
        }

//------------------------------------------------------------------------------

        real Geometry_Interpolator::det_J( const Matrix< DDRMat > & aParamPoint )
        {
            // get tXi and tTau
            Matrix< DDRMat > tXi = aParamPoint( { 0, mNumSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, aParamPoint( mNumSpaceDim ) );

            // get the space jacobian
            Matrix< DDRMat > tdNSpacedXi = this->dNdXi( tXi );
            Matrix< DDRMat > tSpaceJt    = this->space_jacobian( tdNSpacedXi );

            // get the time Jacobian
            Matrix< DDRMat > tdNTimedTau = this->dNdTau( tTau );
            Matrix< DDRMat > tTimeJt     = this->time_jacobian( tdNTimedTau );

            // compute the determinant of the space time Jacobian
            return det( tSpaceJt ) * det( tTimeJt );

        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::time_surf_det_J(       real             & aTimeSurfDetJ,
                                                     const Matrix< DDRMat > & aSideParamPoint,
                                                     const moris_index      & aTimeOrdinal )
        {
            // FIXME time sideset only and LINE, QUAD, HEX only

            // get the space shape function first derivatives wrt the parametric coordinates
            Matrix< DDRMat > tdNSideSpacedXi = this->dNdXi( aSideParamPoint );

            // get the parametric coordinates of the time sideset
            Matrix< DDRMat > tTimeSidesetParamCoords = this->get_time_sideset_param_coords( aTimeOrdinal );

            // evaluation of tangent vectors to the time sideset in the reference space
            Matrix< DDRMat > tTangentVectors = tdNSideSpacedXi * trans( tTimeSidesetParamCoords );

            // get the location of the side integration point in the reference parent element
            Matrix< DDRMat > tParentParamPoint = this->time_surf_val( aSideParamPoint, aTimeOrdinal );

            // get the jacobian
            Matrix< DDRMat > tdNVolSpacedXi = this->dNdXi(  tParentParamPoint( { 0, mNumSpaceDim-1 }, { 0, 0 }) );
            Matrix< DDRMat > tdNVolTimedTau = this->dNdTau( tParentParamPoint( { mNumSpaceDim, mNumSpaceDim }, { 0, 0 }) );
            Matrix< DDRMat > tJt( mNumSpaceDim + mNumTimeDim, mNumSpaceDim + mNumTimeDim, 0.0 );
            tJt( {0,mNumSpaceDim-1}, {0,mNumSpaceDim-1} )
                = this->space_jacobian( tdNVolSpacedXi ).matrix_data();
            tJt( {mNumSpaceDim, mNumSpaceDim}, {mNumSpaceDim, mNumSpaceDim} )
                = this->time_jacobian( tdNVolTimedTau ).matrix_data();

            // evaluation of tangent vectors to the time sideset in the physical space
            Matrix< DDRMat > tRealTangentVectors = trans( tJt ) * trans( tTangentVectors );

            Matrix< DDRMat > tMatrixForDet( mNumSpaceDim + mNumTimeDim, mNumSpaceDim + mNumTimeDim, 1.0 );
            tMatrixForDet({ 1, mNumSpaceDim + mNumTimeDim -1 },{ 0, mNumSpaceDim + mNumTimeDim - 1 })
                = trans( tRealTangentVectors );

            //FIXME abs and /2 because bar1 and not point
            aTimeSurfDetJ =  std::abs( det( tMatrixForDet ) ) / 2.0;

            // FIXME switch on geometry type, cause only valid for line, quad, hex
            // FIXME check that determinant is greater than zero
        }

//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::surf_val( const Matrix< DDRMat > & aSideParamPoint,
                                                           const moris_index      & aSpaceOrdinal )
        {
            // check that there is a side interpolation
            MORIS_ASSERT( mSpaceSideset, "Geometry_Interpolator::surf_val - no side interpolation." );

            // fixme check input size and value

            // get the number of dimensions and bases in space for the side
            uint tSideSpaceDim   = mSideSpaceInterpolation->get_number_of_dimensions();
            uint tSideSpaceBases = mSideSpaceInterpolation->get_number_of_bases();

            // get tXi and tTau
            Matrix< DDRMat > tXi = aSideParamPoint( { 0, tSideSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, aSideParamPoint( tSideSpaceDim ) );

            // evaluate side space interpolation shape functions at aParamPoint
            Matrix< DDRMat > tNSideSpace = mSideSpaceInterpolation->eval_N( tXi );

            // evaluate time interpolation shape functions at aParamPoint
            Matrix< DDRMat > tNTime = mTimeInterpolation->eval_N( tTau );

            // build space time interpolation functions for side
            Matrix< DDRMat > tN = reshape( trans( tNSideSpace ) * tNTime, 1, mNumTimeBases*tSideSpaceBases );

            // get the parametric coordinates of the side in the parent reference element
            Matrix< DDRMat > tSideParentParamCoords = this->get_space_sideset_param_coords( aSpaceOrdinal );

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tParentParamPoint = tN * trans( tSideParentParamCoords );

            return trans( tParentParamPoint );
          }

//------------------------------------------------------------------------------

        Matrix < DDRMat > Geometry_Interpolator::time_surf_val( const Matrix< DDRMat > & aSideParamPoint,
                                                                const moris_index      & aTimeOrdinal )
        {
            // fixme check input size and value

            // get tXi and tTau
            Matrix< DDRMat > tXi = aSideParamPoint( { 0, mNumSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, aSideParamPoint( mNumSpaceDim ) );

            // build space time interpolation functions for side
            Matrix< DDRMat > tN = mSpaceInterpolation->eval_N( tXi );

            // get the parametric coordinates of the side in the parent reference element
            Matrix< DDRMat > tSideParentParamCoords = this->get_time_sideset_param_coords( aTimeOrdinal );

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tParentParamPoint = tN * trans( tSideParentParamCoords );

            //FIXME when computed -1 and 1 are slightly off
            if ( aTimeOrdinal == 0 )
            {
                tParentParamPoint( mNumSpaceDim ) = -1.0;
            }
            else
            {
                tParentParamPoint( mNumSpaceDim ) = 1.0;
            }

            return trans( tParentParamPoint );
        }


//------------------------------------------------------------------------------

        void Geometry_Interpolator::surf_det_J(       real             & aSurfDetJ,
                                                      Matrix< DDRMat > & aNormal,
                                                const Matrix< DDRMat > & aSideParamPoint,
                                                const moris_index      & aSpaceOrdinal )
        {
            // check that there is a side interpolation
            MORIS_ASSERT( mSpaceSideset, "Geometry_Interpolator::surf_val - no side interpolation." );

            // fixme check input size and value

            // get the number of dimensions and bases in space for the side
            uint tSideSpaceDim   = mSideSpaceInterpolation->get_number_of_dimensions();
            uint tSideSpaceBases = mSideSpaceInterpolation->get_number_of_bases();

            // get the number of space time bases
            uint tSideBases = tSideSpaceBases * mNumTimeBases;

            // get tXi and tTau
            Matrix< DDRMat > tXi = aSideParamPoint( { 0, tSideSpaceDim-1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, aSideParamPoint( tSideSpaceDim ) );

            // evaluate side space interpolation shape functions at aParamPoint
            Matrix< DDRMat > tNSideSpace = mSideSpaceInterpolation->eval_N( tXi );

            // evaluate time interpolation shape functions at aParamPoint
            Matrix< DDRMat > tNTime      = mTimeInterpolation->eval_N( tTau );

            // evaluate side space interpolation shape functions first parametric derivatives at aParamPoint
            Matrix< DDRMat > tdNSideSpacedXi = mSideSpaceInterpolation->eval_dNdXi( tXi );
            tdNSideSpacedXi = trans( tdNSideSpacedXi );

            // evaluate time interpolation shape functions first parametric derivatives at aParamPoint
            Matrix< DDRMat > tdNTimedTau     = mTimeInterpolation->eval_dNdXi( tTau );

            // build the space time dNdXi row by row
            Matrix< DDRMat > tdNdXi( tSideSpaceDim, tSideBases );
            for ( uint Ik = 0; Ik < tSideSpaceDim; Ik++ )
            {
                tdNdXi.get_row( Ik ) = reshape( tdNSideSpacedXi.get_column( Ik ) * tNTime , 1, tSideBases );
            }

            // build the space time dNdTau row by row
            Matrix< DDRMat > tdNdTau = reshape( trans( tNSideSpace ) * tdNTimedTau, 1, tSideBases);

            // get the sideset parametric coordinates in the reference volume space
            Matrix< DDRMat > tSideParamCoords = this->get_space_sideset_param_coords( aSpaceOrdinal );

            // evaluation of tangent vectors to the space sideset
            Matrix< DDRMat > tTangentVectors( tSideSpaceDim + mNumTimeDim, mNumSpaceDim + mNumTimeDim );
            tTangentVectors( { 0, tSideSpaceDim-1 }, { 0, mNumSpaceDim + mNumTimeDim-1 } )
                = tdNdXi * trans( tSideParamCoords );
            tTangentVectors( { tSideSpaceDim, tSideSpaceDim }, { 0, mNumSpaceDim + mNumTimeDim-1 } )
                = tdNdTau * trans( tSideParamCoords );
            tTangentVectors = trans( tTangentVectors );

            // initialize the determinant of the Jacobian mapping
            aSurfDetJ = -1.0;
            aNormal.set_size( mNumSpaceDim, 1, 0.0 );

            // get the location of the side integration point in the reference parent element
            Matrix< DDRMat > tParentParamPoint = this->surf_val( aSideParamPoint, aSpaceOrdinal );

            // switch on geometry type
            switch ( mGeometryType )
            {
                case ( mtk::Geometry_Type::LINE ):
                {
                    // get the time jacobian
                    Matrix< DDRMat > tdNTimedTau2 = this->dNdTau( {{tParentParamPoint( 2 )}} );
                    Matrix< DDRMat > tTimeJt = this->time_jacobian( tdNTimedTau2 );

                    aSurfDetJ = tTangentVectors( 1 ) * tTimeJt( 0 );
                    break;
                }

                case ( mtk::Geometry_Type::QUAD ):
                {
                    Matrix< DDRMat > tTRef1 = tTangentVectors.get_column( 0 );
                    Matrix< DDRMat > tTRef2 = tTangentVectors.get_column( 1 );

                    // get the jacobian
                    Matrix< DDRMat > tdNSpacedXi  = this->dNdXi( tParentParamPoint({ 0, 1 },{ 0, 0 }) );
                    Matrix< DDRMat > tdNTimedTau2 = this->dNdTau( {{tParentParamPoint( 2 )}} );

                    Matrix< DDRMat > tJt( 3, 3, 0.0 );
                    tJt( {0,1}, {0,1} ) = this->space_jacobian( tdNSpacedXi ).matrix_data();
                    tJt( {2,2}, {2,2} ) = this->time_jacobian( tdNTimedTau2 ).matrix_data();

                    Matrix< DDRMat > tTReal1 = trans( tJt ) * tTRef1;
                    Matrix< DDRMat > tTReal2 = trans( tJt ) * tTRef2;
                    aSurfDetJ = norm( cross( tTReal1, tTReal2 ) );

//                    Matrix< DDRMat > tMatrixForDet( 3, 3, 1.0 );
//                    tMatrixForDet( {1,1}, {0,2} ) = trans( tTReal1 );
//                    tMatrixForDet( {2,2}, {0,2} ) = trans( tTReal2 );
//                    real aSurfDetJ2 =  std::abs( det( tMatrixForDet ) );
//                    std::cout<<aSurfDetJ2<<std::endl;

                    // FIXME computing the normal
                    aNormal = {{ -tTReal1( 1 ) },
                               {  tTReal1( 0 ) }};
                    aNormal = -1.0 * aNormal / norm( aNormal );
                    //print(tNormal,"tNormal");
                    break;
                }

                case ( mtk::Geometry_Type::HEX ):
                {
                    Matrix< DDRMat > tTRef1 = tTangentVectors.get_column( 0 );
                    Matrix< DDRMat > tTRef2 = tTangentVectors.get_column( 1 );
                    Matrix< DDRMat > tTRef3 = tTangentVectors.get_column( 2 );

                    // get the jacobian
                    Matrix< DDRMat > tdNSpacedXi  = this->dNdXi( tParentParamPoint({ 0, 2 },{ 0, 0 }) );
                    Matrix< DDRMat > tdNTimedTau2 = this->dNdTau( {{tParentParamPoint( 3 )}} );

                    Matrix< DDRMat > tJt( 4, 4, 0.0 );
                    tJt( {0,2}, {0,2} ) = this->space_jacobian( tdNSpacedXi ).matrix_data();
                    tJt( {3,3}, {3,3} ) = this->time_jacobian( tdNTimedTau2 ).matrix_data();

                    Matrix< DDRMat > tTReal1 = trans( tJt ) * tTRef1;
                    Matrix< DDRMat > tTReal2 = trans( tJt ) * tTRef2;
                    Matrix< DDRMat > tTReal3 = trans( tJt ) * tTRef3;

                    Matrix< DDRMat > tMatrixForDet( 4, 4, 1.0 );
                    tMatrixForDet( {1,1}, {0,3} ) = trans( tTReal1 );
                    tMatrixForDet( {2,2}, {0,3} ) = trans( tTReal2 );
                    tMatrixForDet( {3,3}, {0,3} ) = trans( tTReal3 );
                    // FIXME abs?
                    aSurfDetJ = std::abs( det( tMatrixForDet ) );

                    // FIXME normal
                    aNormal = cross( tTReal1({0,2},{0,0}), tTReal2({0,2},{0,0}));
                    aNormal = aNormal / norm( aNormal );
                    //print(aNormal,"aNormal");
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "Geometry_Interpolator::surf_det_J - wrong space dimension");
                    break;
                }

            }
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valx( const Matrix< DDRMat > & aXi )
        {
            //evaluate the field
            return this->NXi( aXi ) * mXHat ;
        }

//------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valt( const Matrix< DDRMat > & aTau )
        {
            //evaluate the field
            return this->NTau( aTau ) * mTHat ;
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives(
                      Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & adNdXi,
                const Matrix< DDRMat > & ad2NdXi2 ) const
        {
            // evaluate transposed of geometry Jacobian
            aJt = this->space_jacobian( adNdXi );

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesSpace( aJt,
                                                  aKt,
                                                  aLt,
                                                  ad2NdXi2,
                                                  mXHat);
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives(
                      Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & adNdTau,
                const Matrix< DDRMat > & ad2NdTau2 ) const
        {
            // evaluate transposed of geometry Jacobian
            aJt = this->time_jacobian( adNdTau );

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesTime( aJt,
                                                 aKt,
                                                 aLt,
                                                 ad2NdTau2,
                                                 mTHat );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat)
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 1, 1 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 3, 3 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
            aLt( 1, 0 ) = std::pow( aJt( 1, 0 ), 2 );
            aLt( 2, 0 ) = aJt( 0 , 0 ) * aJt( 1 , 0 );

            aLt( 0, 1 ) = std::pow( aJt( 0, 1 ), 2 );
            aLt( 1, 1 ) = std::pow( aJt( 1, 1 ), 2 );
            aLt( 2, 1 ) = aJt( 0 , 1 ) * aJt( 1, 1 );

            aLt( 0, 2 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 2 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 2 ) = aJt( 0, 0 )* aJt( 1, 1 ) +  aJt( 0, 1 ) * aJt( 1, 0 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat > & aJt,
                      Matrix< DDRMat > & aKt,
                      Matrix< DDRMat > & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 6, 6 );
            for( uint j=0; j<3; ++j )
            {
                aLt( 0, j ) = std::pow( aJt( 0, j ), 2 );
                aLt( 1, j ) = std::pow( aJt( 1, j ), 2 );
                aLt( 2, j ) = std::pow( aJt( 2, j ), 2 );
                aLt( 3, j ) = aJt( 1 , j ) * aJt( 2 , j );
                aLt( 4, j ) = aJt( 0 , j ) * aJt( 2 , j );
                aLt( 5, j ) = aJt( 0 , j ) * aJt( 1 , j );
            }

            aLt( 0, 3 ) = 2.0 * aJt( 0, 1 ) * aJt( 0, 2 );
            aLt( 1, 3 ) = 2.0 * aJt( 1, 1 ) * aJt( 1, 2 );
            aLt( 2, 3 ) = 2.0 * aJt( 2, 1 ) * aJt( 2, 2 );
            aLt( 3, 3 ) = aJt( 1, 1 ) * aJt( 2, 2 ) + aJt( 2, 1 ) * aJt( 1, 2 );
            aLt( 4, 3 ) = aJt( 0, 1 ) * aJt( 2, 2 ) + aJt( 2, 1 ) * aJt( 0, 2 );
            aLt( 5, 3 ) = aJt( 0, 1 ) * aJt( 1, 2 ) + aJt( 1, 1 ) * aJt( 0, 2 );

            aLt( 0, 4 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 2 );
            aLt( 1, 4 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 2 );
            aLt( 2, 4 ) = 2.0 * aJt( 2, 0 ) * aJt( 2, 2 );
            aLt( 3, 4 ) = aJt( 1, 0 ) * aJt( 2, 2 ) + aJt( 2, 0 ) * aJt( 1, 2 );
            aLt( 4, 4 ) = aJt( 0, 0 ) * aJt( 2, 2 ) + aJt( 2, 0 ) * aJt( 0, 2 );
            aLt( 5, 4 ) = aJt( 0, 0 ) * aJt( 1, 2 ) + aJt( 1, 0 ) * aJt( 0, 2 );

            aLt( 0, 5 ) = 2.0 * aJt( 0, 0 ) * aJt( 0, 1 );
            aLt( 1, 5 ) = 2.0 * aJt( 1, 0 ) * aJt( 1, 1 );
            aLt( 2, 5 ) = 2.0 * aJt( 2, 0 ) * aJt( 2, 1 );
            aLt( 3, 5 ) = aJt( 1, 0 ) * aJt( 2, 1 ) + aJt( 2, 0 ) * aJt( 1, 1 );
            aLt( 4, 5 ) = aJt( 0, 0 ) * aJt( 2, 1 ) + aJt( 2, 0 ) * aJt( 0, 1 );
            aLt( 5, 5 ) = aJt( 0, 0 ) * aJt( 1, 1 ) + aJt( 1, 0 ) * aJt( 0, 1 );
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::set_function_pointers()
        {
            // get number of dimensions and set pointer to function
            // for second space derivative
            switch ( mNumSpaceDim )
            {
                case( 1 ) :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case( 2 ) :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case( 3 ) :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                    break;
                }
            }

            // get number of dimensions and set pointer to function
            // for second time derivative
            switch ( mNumTimeDim )
            {
                case( 1 ) :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case( 2 ) :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case( 3 ) :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::set_function_pointers - unknown number of dimensions. " );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------

        void Geometry_Interpolator::get_auto_side_geometry_type()
        {
            // depending on the parent geometry type
            switch ( mGeometryType )
            {
                // FIXME only QUAD and HEX
                case ( mtk::Geometry_Type::QUAD ):
                {
                    mSideGeometryType = mtk::Geometry_Type::LINE;
                    break;
                }
                case ( mtk::Geometry_Type::HEX ):
                {
                    mSideGeometryType = mtk::Geometry_Type::QUAD;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " Geometry_Interpolator::get_side_geometry_type - undefined geometry type. " );
                    mSideGeometryType = mtk::Geometry_Type::UNDEFINED;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
