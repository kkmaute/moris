#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "fn_dot.hpp"
#include "fn_sum.hpp"
#include "op_div.hpp"
#include "fn_linsolve.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        Geometry_Interpolator::Geometry_Interpolator(
                const Interpolation_Rule & aInterpolationRule,
                const bool                 aSpaceSideset,
                const bool                 aTimeSideset )
        {
            // set bool for side interpolation to true
            mSpaceSideset = aSpaceSideset;

            // set bool for time side interpolation to true
            mTimeSideset = aTimeSideset;

            // create member pointer to space interpolation function
            mSpaceInterpolation = aInterpolationRule.create_space_interpolation_function();

            // create member pointer to time interpolation function
            mTimeInterpolation  = aInterpolationRule.create_time_interpolation_function();

            // number of space bases and dimensions
            mNumSpaceBases    = mSpaceInterpolation->get_number_of_bases();
            mNumSpaceDim      = mSpaceInterpolation->get_number_of_dimensions();
            mNumSpaceParamDim = mSpaceInterpolation->get_number_of_param_dimensions();

            // number of time bases and dimensions
            mNumTimeBases = mTimeInterpolation->get_number_of_bases();
            mNumTimeDim   = mTimeInterpolation->get_number_of_dimensions();

            // set member geometry type
            mGeometryType     = aInterpolationRule.get_geometry_type();
            mTimeGeometryType = aInterpolationRule.get_time_geometry_type();

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
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::reset_eval_flags()
        {
            // reset bool for evaluation
            mNXiEval      = true;
            mdNdXiEval    = true;
            md2NdXi2Eval  = true;
            md3NdXi3Eval  = true;
            mNTauEval     = true;
            mdNdTauEval   = true;
            md2NdTau2Eval = true;
            md3NdTau3Eval = true;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_coeff(
                const Matrix< DDRMat > & aXHat,
                const Matrix< DDRMat > & aTHat )
        {
            // check the space coefficients input size
            // fixme can not check the number of cols for aXHat and aTHat
            MORIS_ASSERT( aXHat.n_rows() == mNumSpaceBases ,
                    " Geometry_Interpolator::set_coeff - Wrong input size (aXHat). ");

            // set the space coefficients
            mXHat = aXHat;

            //check the time coefficients input size
            MORIS_ASSERT( aTHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_coeff - Wrong input size (aTHat). ");

            // set the time coefficients
            mTHat = aTHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space_coeff( const Matrix< DDRMat > & aXHat )
        {
            //check the space coefficients input size
            // fixme can not check the number of cols for aXHat
            MORIS_ASSERT( aXHat.n_rows() == mNumSpaceBases,
                    " Geometry_Interpolator::set_space_coeff - Wrong input size (aXHat). ");

            // set the space coefficients
            mXHat = aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_time_coeff( const Matrix< DDRMat > & aTHat )
        {
            //check the time coefficients input size
            // fixme can not check the number of cols for aTHat
            MORIS_ASSERT( aTHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTHat). ");

            // set the time coefficients
            mTHat = aTHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_param_coeff()
        {
            // default implementation
            // set space and time param coords
            mSpaceInterpolation->get_param_coords( mXiHat );
            mXiHat = trans( mXiHat );

            mTimeInterpolation->get_param_coords( mTauHat );
            mTauHat = trans( mTauHat );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_param_coeff(
                const Matrix< DDRMat > & aXiHat,
                const Matrix< DDRMat > & aTauHat )
        {
            //check the space param coefficients input size
            MORIS_ASSERT( aXiHat.n_rows() == mNumSpaceBases,
                    " Geometry_Interpolator::set_space_param_coeff - Wrong input size (aXiHat). ");

            // set the space coefficients
            mXiHat = aXiHat;

            //check the time param coefficients input size
            MORIS_ASSERT( aTauHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTauHat). ");

            // set the time coefficients
            mTauHat = aTauHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space_param_coeff( const Matrix< DDRMat > & aXiHat )
        {
            //check the space param coefficients input size
            // fixme can not check the number of cols for aXiHat
            MORIS_ASSERT( aXiHat.n_rows() == mNumSpaceBases,
                    " Geometry_Interpolator::set_space_param_coeff - Wrong input size (aXiHat). ");

            // set the space coefficients
            mXiHat = aXiHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_time_param_coeff( const Matrix< DDRMat > & aTauHat )
        {
            //check the time param coefficients input size
            // fixme can not check the number of cols for aTauHat
            MORIS_ASSERT( aTauHat.n_rows() == mNumTimeBases,
                    " Geometry_Interpolator::set_time_coeff - Wrong input size (aTauHat). ");

            // set the time coefficients
            mTauHat = aTauHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space_time( const Matrix< DDRMat > & aParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aParamPoint.n_cols() == 1 ) && ( aParamPoint.n_rows() == mNumSpaceParamDim + mNumTimeDim )),
                    "Geometry_Interpolator::set_space_time - Wrong input size ( aParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumSpaceParamDim + mNumTimeDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aParamPoint( Ik ) <= 1.0 + 1E-12 ) && ( aParamPoint( Ik ) >= -1.0 - 1E-12 ) ),
                        "Geometry_Interpolator::set_space_time - Wrong input value ( aParamPoint ).");
            }

            // set input values
            mXiLocal.matrix_data()  = aParamPoint( { 0, mNumSpaceParamDim-1 }, { 0, 0 } );
            mTauLocal.matrix_data() = aParamPoint( mNumSpaceParamDim );

            // reset bool for evaluation
            this->reset_eval_flags();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_space( const Matrix< DDRMat > & aSpaceParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aSpaceParamPoint.n_cols() == 1 ) && ( aSpaceParamPoint.n_rows() == mNumSpaceParamDim )),
                    "Geometry_Interpolator::set_space - Wrong input size ( aSpaceParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumSpaceParamDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aSpaceParamPoint( Ik ) <= 1.0 + 1E-12 ) && ( aSpaceParamPoint( Ik ) >= -1.0 - 1E-12 ) ),
                        "Geometry_Interpolator::set_space - Wrong input value ( aSpaceParamPoint ).");
            }

            // set input values
            mXiLocal  = aSpaceParamPoint;

            // reset bool for evaluation
            this->reset_eval_flags();
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_time( const Matrix< DDRMat > & aTimeParamPoint )
        {
            // check input size aParamPoint
            MORIS_ASSERT( ( ( aTimeParamPoint.n_cols() == 1 ) && ( aTimeParamPoint.n_rows() == mNumTimeDim )),
                    "Geometry_Interpolator::set_space - Wrong input size ( aTimeParamPoint ).");

            // check input values are between -1 and 1
            // fixme what about TRI and TET
            for ( uint Ik = 0; Ik < mNumTimeDim; Ik++ )
            {
                MORIS_ASSERT( ( ( aTimeParamPoint( Ik ) <= 1.0 + 1E-12 ) && ( aTimeParamPoint( Ik ) >= -1.0 - 1E-12 ) ),
                        "Geometry_Interpolator::set_time - Wrong input value ( aTimeParamPoint ).");
            }

            // set input values
            mTauLocal  = aTimeParamPoint;

            // reset bool for evaluation
            this->reset_eval_flags();
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::get_time_step()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel()>0, "Geometry_Interpolator::get_time_step - mTHat is not set." );

            // compute time increment deltat
            return mTHat.max() - mTHat.min();
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::extract_space_side_space_param_coeff(
                moris_index              aSpaceOrdinal,
                mtk::Interpolation_Order aInterpolationOrder )
        {
            // get the interpolation order of the IP cell
            mtk::Interpolation_Order tIPInterpolationOrder = mSpaceInterpolation->get_interpolation_order();

            Matrix< DDRMat > tXiForExtraction;
            uint tNumSpaceBases;

            if( tIPInterpolationOrder == aInterpolationOrder )
            {
                // get the parametric coordinated of the parent for extraction
                mSpaceInterpolation->get_param_coords( tXiForExtraction );
                tXiForExtraction = trans( tXiForExtraction );

                // get the number of bases
                tNumSpaceBases = mNumSpaceBases;
            }
            else
            {
                // create an interpolation rule for the iG cell
                Interpolation_Rule tIGInterpolationRule( mGeometryType,
                        Interpolation_Type::LAGRANGE,
                        aInterpolationOrder,
                        Interpolation_Type::LAGRANGE,
                        mTimeInterpolation->get_interpolation_order() );

                // create an Interpolation function for the IG cell
                Interpolation_Function_Base * tIGSpaceInterpolation = tIGInterpolationRule.create_space_interpolation_function();

                // get the parametric coordinated of the parent for extraction
                tIGSpaceInterpolation->get_param_coords( tXiForExtraction );
                tXiForExtraction = trans( tXiForExtraction );

                // get the number of bases
                tNumSpaceBases = tIGSpaceInterpolation->get_number_of_bases();

                // clean up
                delete tIGSpaceInterpolation;
            }

            // get the vertices per side map
            moris::Cell< moris::Cell < moris_index > >tVerticesPerFace =
                    this->get_face_vertices_ordinals( tNumSpaceBases );

            // check if the space ordinal is appropriate
            MORIS_ASSERT( aSpaceOrdinal < static_cast< moris_index > ( tVerticesPerFace.size() ) ,
                    "Geometry_Interpolator::extract_space_side_space_param_coeff - wrong ordinal" );

            // get space side vertices ordinals
            moris::Cell< moris::moris_index > tVerticesOrdinals = tVerticesPerFace( aSpaceOrdinal );

            // initialize the time side set parametric coordinates matrix
            uint tNumOfVertices = tVerticesOrdinals.size();

            // init the matrix with the face param coords
            Matrix< DDRMat > tSideXiHat( tNumOfVertices, mNumSpaceParamDim );

            // loop over the vertices of the face
            for( uint i = 0; i < tNumOfVertices; i++ )
            {
                // get the treated vertex
                moris_index tTreatedVertex = tVerticesOrdinals( i );

                // get the vertex parametric coordinates for node tTreatedVertex
                tSideXiHat.get_row( i ) = tXiForExtraction.get_row( tTreatedVertex );
            }
            return tSideXiHat;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::extract_space_param_coeff(
                mtk::Interpolation_Order aInterpolationOrder )
        {
            // get the interpolation order of the IP cell
            mtk::Interpolation_Order tIPInterpolationOrder =
                    mSpaceInterpolation->get_interpolation_order();

            Matrix< DDRMat > tXi;
            if( tIPInterpolationOrder == aInterpolationOrder )
            {
                // get the parametric coordinated of the parent for extraction
                mSpaceInterpolation->get_param_coords( tXi );
                tXi = trans( tXi );
            }
            else
            {
                // create an interpolation rule for the iG cell
                Interpolation_Rule tIGInterpolationRule(
                        mGeometryType,
                        Interpolation_Type::LAGRANGE,
                        aInterpolationOrder,
                        Interpolation_Type::LAGRANGE,
                        mTimeInterpolation->get_interpolation_order() );

                // create an Interpolation function for the IG cell
                Interpolation_Function_Base * tIGSpaceInterpolation =
                        tIGInterpolationRule.create_space_interpolation_function();

                // get the parametric coordinated of the parent for extraction
                tIGSpaceInterpolation->get_param_coords( tXi );
                tXi = trans( tXi );

                // clean up
                delete tIGSpaceInterpolation;
            }
            return tXi;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::NXi()
        {
            // if shape functions need to be evaluated
            if( mNXiEval )
            {
                // evaluate the shape functions
                this->eval_NXi();

                // set bool for evaluation
                mNXiEval = false;
            }

            // return member value
            return mNXi;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_NXi()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_NXi - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_N( mXiLocal, mNXi );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::NTau()
        {
            // if shape functions need to be evaluated
            if( mNTauEval )
            {
                // evaluate the shape functions
                this->eval_NTau();

                // set bool for evaluation
                mNTauEval = false;
            }

            // return member value
            return mNTau;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_NTau()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_NTau - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_N( mTauLocal, mNTau );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::dNdXi()
        {
            // if shape functions need to be evaluated
            if( mdNdXiEval )
            {
                // evaluate the shape functions 1st derivative
                this->eval_dNdXi();

                // set bool for evaluation
                mdNdXiEval = false;
            }

            // return member value
            return mdNdXi;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_dNdXi()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_dNdXi - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_dNdXi( mXiLocal, mdNdXi );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::dNdTau()
        {
            // if shape functions need to be evaluated
            if( mdNdTauEval )
            {
                // evaluate the shape functions 1st derivative
                this->eval_dNdTau();

                // set bool for evaluation
                mdNdTauEval = false;
            }

            // return member value
            return mdNdTau;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_dNdTau()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_dNdTau - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_dNdXi( mTauLocal, mdNdTau );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::d2NdXi2()
        {
            // if shape functions need to be evaluated
            if( md2NdXi2Eval )
            {
                // evaluate the shape functions 2nd derivative
                this->eval_d2NdXi2();

                // set bool for evaluation
                md2NdXi2Eval = false;
            }

            // return member value
            return md2NdXi2;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_d2NdXi2()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d2NdXi2 - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_d2NdXi2( mXiLocal, md2NdXi2 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::d3NdXi3()
        {
            // if shape functions need to be evaluated
            if( md3NdXi3Eval )
            {
                // evaluate the shape functions 3rd derivative
                this->eval_d3NdXi3();

                // set bool for evaluation
                md3NdXi3Eval = false;
            }

            // return member value
            return md3NdXi3;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_d3NdXi3()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d3NdXi3 - mXiLocal is not set." );

            // pass data through interpolation function
            mSpaceInterpolation->eval_d3NdXi3( mXiLocal, md3NdXi3 );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Geometry_Interpolator::d2NdTau2()
        {
            // if shape functions need to be evaluated
            if( md2NdTau2Eval )
            {
                // evaluate the shape functions 2nd derivative
                this->eval_d2NdTau2();

                // set bool for evaluation
                md2NdTau2Eval = false;
            }

            // return member value
            return md2NdTau2;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_d2NdTau2()
        {
            // check that mXiLocal is set
            MORIS_ASSERT( mTauLocal.numel() > 0,
                    "Geometry_Interpolator::eval_d2NdTau2 - mTauLocal is not set." );

            // pass data through interpolation function
            mTimeInterpolation->eval_d2NdXi2( mTauLocal, md2NdTau2 );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian( Matrix< DDRMat > & aJt )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel()>0,
                    "Geometry_Interpolator::space_jacobian - mXHat is not set." );

            // compute the Jacobian
            aJt = this->dNdXi() * mXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::second_space_jacobian( Matrix< DDRMat > & aJ2bt )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel()>0,
                    "Geometry_Interpolator::second_space_jacobian - mXHat is not set." );

            // compute the second order Jacobian
            aJ2bt = this->d2NdXi2() * mXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::third_space_jacobian( Matrix< DDRMat > & aJ3ct )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel()>0,
                    "Geometry_Interpolator::third_space_jacobian - mXHat is not set." );

            // compute the third order Jacobian
            aJ3ct = this->d3NdXi3() * mXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::time_jacobian( Matrix< DDRMat > & aJt )
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel()>0,
                    "Geometry_Interpolator::time_jacobian - mTHat is not set." );

            // compute the Jacobian
            aJt = this->dNdTau() * mTHat;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::det_J()
        {
            // get the space jacobian
            Matrix< DDRMat > tSpaceJt;
            this->space_jacobian( tSpaceJt );

            // get the time Jacobian
            Matrix< DDRMat > tTimeJt;
            this->time_jacobian( tTimeJt );

            // compute detJ for space
            real detJSpace = ( this->*mSpaceDetJFunc )( tSpaceJt );

            // compute detJ for time
            real detJTime  = ( this->*mTimeDetJFunc )( tTimeJt );

            // compute the determinant of the space time Jacobian
            return detJSpace * detJTime;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_side_line(
                Matrix< DDRMat > & aSpaceJt )
        {
            return norm( aSpaceJt );
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_side_tri(
                Matrix< DDRMat > & aSpaceJt )
        {
            return norm( cross(
                    aSpaceJt.get_row( 0 ) - aSpaceJt.get_row( 2 ),
                    aSpaceJt.get_row( 1 ) - aSpaceJt.get_row( 2 ) ) ) / 2.0;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_side_quad(
                Matrix< DDRMat > & aSpaceJt )
        {
            return norm( cross( aSpaceJt.get_row( 0 ), aSpaceJt.get_row( 1 ) ) );
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_bulk_line_quad_hex(
                Matrix< DDRMat > & aSpaceJt )
        {
            return det( aSpaceJt );
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_bulk_tri(
                Matrix< DDRMat > & aSpaceJt )
        {
            Matrix< DDRMat > tSpaceJt2( mNumSpaceParamDim, mNumSpaceParamDim, 1.0 );

            tSpaceJt2( { 1, mNumSpaceParamDim - 1 },{ 0, mNumSpaceParamDim - 1 } ) =
                    trans( aSpaceJt );

            return det( tSpaceJt2 ) / 2.0;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_space_detJ_bulk_tet(
                Matrix< DDRMat > & aSpaceJt )
        {
            Matrix< DDRMat > tSpaceJt2( mNumSpaceParamDim, mNumSpaceParamDim, 1.0 );

            tSpaceJt2( { 1, mNumSpaceParamDim - 1 }, { 0, mNumSpaceParamDim - 1 } ) =
                    trans( aSpaceJt );

            return det( tSpaceJt2 ) / 6.0;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_time_detJ_side(
                Matrix< DDRMat > & aTimeJt )
        {
            return 1.0;
        }

        //------------------------------------------------------------------------------

        real Geometry_Interpolator::eval_time_detJ_bulk(
                Matrix< DDRMat > & aTimeJt )
        {
            return det( aTimeJt );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::get_normal( Matrix< DDRMat > & aNormal )
        {
            // check that there is a side interpolation
            MORIS_ASSERT( mSpaceSideset,
                    "Geometry_Interpolator::normal - not a side." );

            // check that mXiLocal is set
            MORIS_ASSERT( mXiLocal.numel() > 0,
                    "Geometry_Interpolator::normal - mXiLocal is not set." );

            // evaluate side space interpolation shape functions first parametric derivatives at aParamPoint
            Matrix< DDRMat > tdNSpacedXi;
            mSpaceInterpolation->eval_dNdXi( mXiLocal, tdNSpacedXi );

            // evaluation of tangent vectors to the space side in the physical space
            Matrix< DDRMat > tRealTangents = trans( tdNSpacedXi * mXHat );

            // computing the normal from the real tangent vectors
            ( this->*mNormalFunc )( tRealTangents, aNormal );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_normal_side_line(
                Matrix< DDRMat > & aRealTangents,
                Matrix< DDRMat > & aNormal )
        {
            aNormal = {
                    {   aRealTangents( 1 ) },
                    { - aRealTangents( 0 ) }};

            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_normal_side_quad(
                Matrix< DDRMat > & aRealTangents,
                Matrix< DDRMat > & aNormal )
        {
            aNormal = cross( aRealTangents.get_column( 0 ), aRealTangents.get_column( 1 ) );
            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_normal_side_tri(
                Matrix< DDRMat > & aRealTangents,
                Matrix< DDRMat > & aNormal )
        {
            aNormal = cross(
                    aRealTangents.get_column( 0 ) - aRealTangents.get_column( 2 ),
                    aRealTangents.get_column( 1 ) - aRealTangents.get_column( 2 ) );
            aNormal = aNormal / norm( aNormal );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valx()
        {
            // check that mTHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::valx - mXHat is not set." );

            //evaluate the field
            return this->NXi() * mXHat ;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Geometry_Interpolator::valt()
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel()>0,
                    "Geometry_Interpolator::valt - mTHat is not set." );

            //evaluate the field
            return this->NTau() * mTHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::map_integration_point(
                Matrix< DDRMat > & aGlobalParamPoint )
        {
            // check that mXiHat and mTauHat are set
            MORIS_ASSERT( mXiHat.numel() > 0,
                    "Geometry_Interpolator::map_integration_point - mXiHat is not set." );
            MORIS_ASSERT( mTauHat.numel() > 0,
                    "Geometry_Interpolator::map_integration_point - mTauHat is not set." );

            // evaluate the coords of the mapped param point
            uint tNumSpaceCoords = mXiHat.n_cols();

            aGlobalParamPoint.set_size( tNumSpaceCoords + 1, 1 );

            aGlobalParamPoint( { 0, tNumSpaceCoords - 1 }, { 0, 0 } ) =
                    trans( this->NXi()  * mXiHat );

            aGlobalParamPoint( { tNumSpaceCoords, tNumSpaceCoords }, { 0, 0 } ) =
                    trans( this->NTau() * mTauHat ) ;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::update_local_coordinates(
                Matrix< DDRMat > & aPhysCoordinates,
                Matrix< DDRMat > & aParamCoordinates )
        {
            // set max iteration
            // std::cout<<"Geometry_Interpolator::update_local_coordinates - Force 10 iterations"<<std::endl;
            uint tMaxIter = 10;

            // init previous param coords
            Matrix< DDRMat > tPreviousParamCoordinates(
                    aParamCoordinates.n_rows(),
                    aParamCoordinates.n_cols(),
                    10.0 );

            // Newton loop
            for( uint iIter = 0; iIter <= tMaxIter; iIter++ )
            {
                // break update
                real tError = norm( tPreviousParamCoordinates - aParamCoordinates );
                if( tError < 1e-12 )
                {
                    break;
                }

                // compute NXi
                Matrix< DDRMat > tNXi;
                mSpaceInterpolation->eval_N( aParamCoordinates, tNXi );

                // compute residual
                Matrix< DDRMat > tR = trans( aPhysCoordinates - tNXi * mXHat );

                // compute dNdXI
                Matrix< DDRMat > tdNdXi;
                mSpaceInterpolation->eval_dNdXi( aParamCoordinates, tdNdXi );

                // compute jacobian
                Matrix< DDRMat > tJ = trans( tdNdXi * mXHat );

                // solve
                aParamCoordinates.matrix_data() += trans( solve( tJ, tR ) );

                // update previous coords
                tPreviousParamCoordinates = aParamCoordinates;

                if( iIter == tMaxIter )
                {
                    MORIS_LOG_INFO("Geometry_Interpolator::update_local_coordinates - No convergence.");
                }
            }
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >       & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & adNdXi,
                const Matrix< DDRMat > & ad2NdXi2 )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::space_jacobian_and_matrices_for_second_derivatives - mXHat is not set." );

            // evaluate transposed of geometry Jacobian
            this->space_jacobian( aJt );

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesSpace(
                    aJt, // contains first geometric derivs
                    aKt,
                    aLt,
                    ad2NdXi2,
                    mXHat );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::space_jacobian_and_matrices_for_third_derivatives(
                Matrix< DDRMat >       & aJt,   // contains first geometric derivs
                Matrix< DDRMat >       & aJ2bt, // contains second geometric derivs = second help matrix for 2nd field derivs
                Matrix< DDRMat >       & aJ3at, // first help matrix for 3rd field derivs
                Matrix< DDRMat >       & aJ3bt, // second help matrix for 3rd field derivs
                Matrix< DDRMat >       & aJ3ct, // third help matrix for 3rd field derivs
                const Matrix< DDRMat > & adNdXi,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & ad3NdXi3 )
        {
            // check that mXHat is set
            MORIS_ASSERT( mXHat.numel() > 0,
                    "Geometry_Interpolator::space_jacobian_and_matrices_for_third_derivatives - mXHat is not set." );

            // evaluate geometry Jacobians
            this->space_jacobian( aJt );
            this->second_space_jacobian( aJ2bt );

            // call calculator for second derivatives
            this->mThirdDerivativeMatricesSpace(
                    aJt,
                    aJ2bt,
                    aJ3at,
                    aJ3bt,
                    aJ3ct,
                    ad3NdXi3,
                    mXHat );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives(
                Matrix< DDRMat >       & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & adNdTau,
                const Matrix< DDRMat > & ad2NdTau2 )
        {
            // check that mTHat is set
            MORIS_ASSERT( mTHat.numel() > 0,
                    "Geometry_Interpolator::time_jacobian_and_matrices_for_second_derivatives - mTHat is not set." );

            // evaluate transposed of geometry Jacobian
            this->time_jacobian( aJt );

            // call calculator for second derivatives
            this->mSecondDerivativeMatricesTime(
                    aJt,
                    aKt,
                    aLt,
                    ad2NdTau2,
                    mTHat );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_1d(
                const Matrix< DDRMat > & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 1, 1 );
            aLt( 0, 0 ) = std::pow( aJt( 0, 0 ), 2 );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_third_derivative_1d(
                const Matrix< DDRMat > & aJt,
                const Matrix< DDRMat > & aJ2bt,
                Matrix< DDRMat >       & aJ3at,
                Matrix< DDRMat >       & aJ3bt,
                Matrix< DDRMat >       & aJ3ct,
                const Matrix< DDRMat > & ad3NdXi3,
                const Matrix< DDRMat > & aXHat )
        {
            // first help matrix
            aJ3at.set_size( 1, 1 );
            aJ3at( 0, 0 ) = std::pow( aJt( 0, 0 ), 3 );

            // second help matrix
            aJ3bt.set_size( 1, 1 );
            aJ3bt( 0, 0 ) = 3 * aJ2bt( 0, 0 ) * aJt( 0, 0 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_2d(
                const Matrix< DDRMat > & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
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
            aLt( 2, 2 ) = aJt( 0, 0 ) * aJt( 1, 1 ) + aJt( 0, 1 ) * aJt( 1, 0 );
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_third_derivative_2d(
                const Matrix< DDRMat > & aJt,
                const Matrix< DDRMat > & aJ2bt,
                Matrix< DDRMat >       & aJ3at,
                Matrix< DDRMat >       & aJ3bt,
                Matrix< DDRMat >       & aJ3ct,
                const Matrix< DDRMat > & ad3NdXi3,
                const Matrix< DDRMat > & aXHat )
        {
            // first help matrix
            aJ3at.set_size( 4, 4 );

            /* matrix structured into 4 parts
             *  _____________     ________
             *  |(1)* |(2)* |     | ,xxx |
             *  |_*_*_|_*_*_|  *  | ,yyy |
             *  |(3)* |(4)* |     | ,xxy |
             *  |_*_*_|_*_*_|     |_,xyy_|
             */

            // Block (1) ------------------------------------------------
            for( uint j = 0; j < 2; ++j )
            {
                aJ3at( 0, j ) = std::pow( aJt( 0, j ), 3 );
                aJ3at( 1, j ) = std::pow( aJt( 1, j ), 3 );
            }

            // Block (2) ------------------------------------------------
            aJ3at( 0, 2 ) = 3 * std::pow( aJt( 0, 0 ), 2 ) * aJt( 0, 1 );
            aJ3at( 1, 2 ) = 3 * std::pow( aJt( 1, 0 ), 2 ) * aJt( 1, 1 );

            aJ3at( 0, 3 ) = 3 * std::pow( aJt( 0, 1 ), 2 ) * aJt( 0, 0 );
            aJ3at( 1, 3 ) = 3 * std::pow( aJt( 1, 1 ), 2 ) * aJt( 1, 0 );

            // Block (3) ------------------------------------------------
            for( uint j = 0; j < 2; ++j )
            {
                aJ3at( 2, j ) = std::pow( aJt( 0, j ), 2 ) * aJt( 1, j );
                aJ3at( 3, j ) = std::pow( aJt( 1, j ), 2 ) * aJt( 0, j );
            }

            // Block (4) ------------------------------------------------
            aJ3at( 2, 2 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 1, 1 )  +  2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 1 );
            aJ3at( 3, 2 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 0, 1 )  +  2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 1 );

            aJ3at( 2, 3 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 1, 0 )  +  2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 0 );
            aJ3at( 3, 3 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 0, 0 )  +  2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 0 );

            // second help matrix
            aJ3bt.set_size( 4, 3 );

            /* matrix structured into 4 parts
             *  ___________     _______
             *  |(1)* |(2)|     | ,xx |
             *  |_*_*_|_*_|  *  | ,yy |
             *  |(3)* |(4)|     |_,xy_|
             *  |_*_*_|_*_|
             */

            // Block (1) ------------------------------------------------
            for( uint j=0; j<2; ++j )
            {
                aJ3bt( 0, j ) = 3 * aJ2bt( 0, j ) * aJt( 0, j );
                aJ3bt( 1, j ) = 3 * aJ2bt( 1, j ) * aJt( 1, j );
            }

            // Block (2) ------------------------------------------------
            aJ3bt( 0, 2 ) =   3 * aJ2bt( 0, 0 ) * aJt( 0, 1 ) + 3 * aJ2bt( 0, 1 ) * aJt( 0, 0 );
            aJ3bt( 1, 2 ) =   3 * aJ2bt( 1, 0 ) * aJt( 1, 1 ) + 3 * aJ2bt( 1, 1 ) * aJt( 1, 0 );

            // Block (3) ------------------------------------------------
            for( uint j=0; j<2; ++j )
            {
                aJ3bt( 2, j ) = 2 * aJ2bt( 2, j ) * aJt( 0, j )  +  aJ2bt( 0, j ) * aJt( 1, j );
                aJ3bt( 3, j ) = 2 * aJ2bt( 2, j ) * aJt( 1, j )  +  aJ2bt( 1, j ) * aJt( 0, j );
            }

            // Block (4) ------------------------------------------------
            aJ3bt( 2, 2 ) =
                    2 * aJ2bt( 2, 0 ) * aJt( 0, 1 )  +  2 * aJ2bt( 2, 1 ) * aJt( 0, 0 ) +
                    1 * aJ2bt( 0, 1 ) * aJt( 1, 0 )  +  1 * aJ2bt( 0, 0 ) * aJt( 1, 1 );
            aJ3bt( 3, 2 ) =
                    2 * aJ2bt( 2, 0 ) * aJt( 1, 1 )  +  2 * aJ2bt( 2, 1 ) * aJt( 1, 0 ) +
                    1 * aJ2bt( 1, 1 ) * aJt( 0, 0 )  +  1 * aJ2bt( 1, 0 ) * aJt( 0, 1 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::eval_matrices_for_second_derivative_3d(
                const Matrix< DDRMat > & aJt,
                Matrix< DDRMat >       & aKt,
                Matrix< DDRMat >       & aLt,
                const Matrix< DDRMat > & ad2NdXi2,
                const Matrix< DDRMat > & aXHat )
        {
            // help matrix K
            aKt = ad2NdXi2 * aXHat;

            // help matrix L
            aLt.set_size( 6, 6 );
            for( uint j = 0; j < 3; ++j )
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

        void Geometry_Interpolator::eval_matrices_for_third_derivative_3d(
                const Matrix< DDRMat > & aJt,
                const Matrix< DDRMat > & aJ2bt,
                Matrix< DDRMat >       & aJ3at,
                Matrix< DDRMat >       & aJ3bt,
                Matrix< DDRMat >       & aJ3ct,
                const Matrix< DDRMat > & ad3NdXi3,
                const Matrix< DDRMat > & aXHat )
        {
            // first help matrix
            aJ3at.set_size( 10, 10 );

            /* matrix structured into 9 parts
             *  ___________________________     ________
             *  | * * * | * * * * * * | * |     | ,xxx |
             *  | *(1)* | * *(2)* * * |(3)|     | ,yyy |
             *  |_*_*_*_|_*_*_*_*_*_*_|_*_|     | ,zzz |
             *  | * * * | * * * * * * | * |     | ,xxy |
             *  | * * * | * * * * * * | * |     | ,xxz |
             *  | *(4)* | * *(5)* * * |(6)|  *  | ,xyy |
             *  | * * * | * * * * * * | * |     | ,yyz |
             *  | * * * | * * * * * * | * |     | ,xzz |
             *  |_*_*_*_|_*_*_*_*_*_*_|_*_|     | ,yzz |
             *  |_*(7)*_|_*_*(8)*_*_*_|(9)|     |_,xyz_|
             */

            // Block (1) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3at( 0, j ) = std::pow( aJt( 0, j ), 3 );
                aJ3at( 1, j ) = std::pow( aJt( 1, j ), 3 );
                aJ3at( 2, j ) = std::pow( aJt( 2, j ), 3 );
            }

            // Block (2) ------------------------------------------------
            aJ3at( 0, 3 ) = 3 * std::pow( aJt( 0, 0 ), 2 ) * aJt( 0, 1 );
            aJ3at( 1, 3 ) = 3 * std::pow( aJt( 1, 0 ), 2 ) * aJt( 1, 1 );
            aJ3at( 2, 3 ) = 3 * std::pow( aJt( 2, 0 ), 2 ) * aJt( 2, 1 );

            aJ3at( 0, 4 ) = 3 * std::pow( aJt( 0, 0 ), 2 ) * aJt( 0, 2 );
            aJ3at( 1, 4 ) = 3 * std::pow( aJt( 1, 0 ), 2 ) * aJt( 1, 2 );
            aJ3at( 2, 4 ) = 3 * std::pow( aJt( 2, 0 ), 2 ) * aJt( 2, 2 );

            aJ3at( 0, 5 ) = 3 * std::pow( aJt( 0, 1 ), 2 ) * aJt( 0, 0 );
            aJ3at( 1, 5 ) = 3 * std::pow( aJt( 1, 1 ), 2 ) * aJt( 1, 0 );
            aJ3at( 2, 5 ) = 3 * std::pow( aJt( 2, 1 ), 2 ) * aJt( 2, 0 );

            aJ3at( 0, 6 ) = 3 * std::pow( aJt( 0, 1 ), 2 ) * aJt( 0, 2 );
            aJ3at( 1, 6 ) = 3 * std::pow( aJt( 1, 1 ), 2 ) * aJt( 1, 2 );
            aJ3at( 2, 6 ) = 3 * std::pow( aJt( 2, 1 ), 2 ) * aJt( 2, 2 );

            aJ3at( 0, 7 ) = 3 * std::pow( aJt( 0, 2 ), 2 ) * aJt( 0, 0 );
            aJ3at( 1, 7 ) = 3 * std::pow( aJt( 1, 2 ), 2 ) * aJt( 1, 0 );
            aJ3at( 2, 7 ) = 3 * std::pow( aJt( 2, 2 ), 2 ) * aJt( 2, 0 );

            aJ3at( 0, 8 ) = 3 * std::pow( aJt( 0, 2 ), 2 ) * aJt( 0, 1 );
            aJ3at( 1, 8 ) = 3 * std::pow( aJt( 1, 2 ), 2 ) * aJt( 1, 1 );
            aJ3at( 2, 8 ) = 3 * std::pow( aJt( 2, 2 ), 2 ) * aJt( 2, 1 );


            // Block (3) ------------------------------------------------
            aJ3at( 0, 9 ) = 6 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 1, 9 ) = 6 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 2, 9 ) = 6 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );

            // Block (4) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3at( 3, j ) = std::pow( aJt( 0, j ) , 2 ) * aJt( 1, j );
                aJ3at( 4, j ) = std::pow( aJt( 0, j ) , 2 ) * aJt( 2, j );
                aJ3at( 5, j ) = std::pow( aJt( 1, j ) , 2 ) * aJt( 0, j );
                aJ3at( 6, j ) = std::pow( aJt( 1, j ) , 2 ) * aJt( 2, j );
                aJ3at( 7, j ) = std::pow( aJt( 2, j ) , 2 ) * aJt( 0, j );
                aJ3at( 8, j ) = std::pow( aJt( 2, j ) , 2 ) * aJt( 1, j );
            }

            // Block (5) ------------------------------------------------
            aJ3at( 3, 3 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 1 );
            aJ3at( 4, 3 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 0, 0 ) * aJt( 2, 0 ) * aJt( 0, 1 );
            aJ3at( 5, 3 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 1 );
            aJ3at( 6, 3 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 1, 0 ) * aJt( 2, 0 ) * aJt( 1, 1 );
            aJ3at( 7, 3 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 2, 0 ) * aJt( 0, 0 ) * aJt( 2, 1 );
            aJ3at( 8, 3 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 2, 0 ) * aJt( 1, 0 ) * aJt( 2, 1 );

            aJ3at( 3, 4 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 0, 2 );
            aJ3at( 4, 4 ) = std::pow( aJt( 0, 0 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 0, 0 ) * aJt( 2, 0 ) * aJt( 0, 2 );
            aJ3at( 5, 4 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 1, 0 ) * aJt( 0, 0 ) * aJt( 1, 2 );
            aJ3at( 6, 4 ) = std::pow( aJt( 1, 0 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 1, 0 ) * aJt( 2, 0 ) * aJt( 1, 2 );
            aJ3at( 7, 4 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 2, 0 ) * aJt( 0, 0 ) * aJt( 2, 2 );
            aJ3at( 8, 4 ) = std::pow( aJt( 2, 0 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 2, 0 ) * aJt( 1, 0 ) * aJt( 2, 2 );

            aJ3at( 3, 5 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 0 );
            aJ3at( 4, 5 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 0, 1 ) * aJt( 2, 1 ) * aJt( 0, 0 );
            aJ3at( 5, 5 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 0 );
            aJ3at( 6, 5 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 1, 1 ) * aJt( 2, 1 ) * aJt( 1, 0 );
            aJ3at( 7, 5 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 2, 1 ) * aJt( 0, 1 ) * aJt( 2, 0 );
            aJ3at( 8, 5 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 2, 1 ) * aJt( 1, 1 ) * aJt( 2, 0 );

            aJ3at( 3, 6 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 0, 2 );
            aJ3at( 4, 6 ) = std::pow( aJt( 0, 1 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 0, 1 ) * aJt( 2, 1 ) * aJt( 0, 2 );
            aJ3at( 5, 6 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 1, 1 ) * aJt( 0, 1 ) * aJt( 1, 2 );
            aJ3at( 6, 6 ) = std::pow( aJt( 1, 1 ), 2 ) * aJt( 2, 2 ) + 2 * aJt( 1, 1 ) * aJt( 2, 1 ) * aJt( 1, 2 );
            aJ3at( 7, 6 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 0, 2 ) + 2 * aJt( 2, 1 ) * aJt( 0, 1 ) * aJt( 2, 2 );
            aJ3at( 8, 6 ) = std::pow( aJt( 2, 1 ), 2 ) * aJt( 1, 2 ) + 2 * aJt( 2, 1 ) * aJt( 1, 1 ) * aJt( 2, 2 );

            aJ3at( 3, 7 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 0, 0 );
            aJ3at( 4, 7 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 0, 2 ) * aJt( 2, 2 ) * aJt( 0, 0 );
            aJ3at( 5, 7 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 1, 2 ) * aJt( 0, 2 ) * aJt( 1, 0 );
            aJ3at( 6, 7 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 2, 0 ) + 2 * aJt( 1, 2 ) * aJt( 2, 2 ) * aJt( 1, 0 );
            aJ3at( 7, 7 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 0, 0 ) + 2 * aJt( 2, 2 ) * aJt( 0, 2 ) * aJt( 2, 0 );
            aJ3at( 8, 7 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 1, 0 ) + 2 * aJt( 2, 2 ) * aJt( 1, 2 ) * aJt( 2, 0 );

            aJ3at( 3, 8 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 0, 1 );
            aJ3at( 4, 8 ) = std::pow( aJt( 0, 2 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 0, 2 ) * aJt( 2, 2 ) * aJt( 0, 1 );
            aJ3at( 5, 8 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 1, 2 ) * aJt( 0, 2 ) * aJt( 1, 1 );
            aJ3at( 6, 8 ) = std::pow( aJt( 1, 2 ), 2 ) * aJt( 2, 1 ) + 2 * aJt( 1, 2 ) * aJt( 2, 2 ) * aJt( 1, 1 );
            aJ3at( 7, 8 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 0, 1 ) + 2 * aJt( 2, 2 ) * aJt( 0, 2 ) * aJt( 2, 1 );
            aJ3at( 8, 8 ) = std::pow( aJt( 2, 2 ), 2 ) * aJt( 1, 1 ) + 2 * aJt( 2, 2 ) * aJt( 1, 2 ) * aJt( 2, 1 );


            // Block (6) ------------------------------------------------
            aJ3at( 3, 9 ) =
                    2 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 4, 9 ) =
                    2 * aJt( 0, 0 ) * aJt( 0, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 2, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 0, 1 ) * aJt( 0, 2 );
            aJ3at( 5, 9 ) =
                    2 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 0, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 6, 9 ) =
                    2 * aJt( 1, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 2, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 1, 1 ) * aJt( 1, 2 );
            aJ3at( 7, 9 ) =
                    2 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 0, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 0, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 0, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );
            aJ3at( 8, 9 ) =
                    2 * aJt( 2, 0 ) * aJt( 2, 1 ) * aJt( 1, 2 ) +
                    2 * aJt( 2, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    2 * aJt( 1, 0 ) * aJt( 2, 1 ) * aJt( 2, 2 );

            // Block (7) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3at( 9, j ) = aJt( 0, j ) * aJt( 1, j ) * aJt( 2, j );
            }

            // Block (8) ------------------------------------------------
            aJ3at( 9, 3 ) =
                    aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 0 ) ;

            aJ3at( 9, 4 ) =
                    aJt( 0, 0 ) * aJt( 1, 0 ) * aJt( 2, 2 ) +
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 0 );

            aJ3at( 9, 5 ) =
                    aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 1 );

            aJ3at( 9, 6 ) =
                    aJt( 0, 1 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 1 );

            aJ3at( 9, 7 ) =
                    aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 2 ) +
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 2 );

            aJ3at( 9, 8 ) =
                    aJt( 0, 2 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 2 );

            // Block (9) ------------------------------------------------
            aJ3at( 9, 9 ) =
                    aJt( 0, 0 ) * aJt( 1, 1 ) * aJt( 2, 2 ) +
                    aJt( 0, 2 ) * aJt( 1, 1 ) * aJt( 2, 0 ) +
                    aJt( 0, 1 ) * aJt( 1, 2 ) * aJt( 2, 0 ) +
                    aJt( 0, 0 ) * aJt( 1, 2 ) * aJt( 2, 1 ) +
                    aJt( 0, 2 ) * aJt( 1, 0 ) * aJt( 2, 1 ) +
                    aJt( 0, 1 ) * aJt( 1, 0 ) * aJt( 2, 2 );

            // second help matrix
            aJ3bt.set_size( 10, 6 );

            /* matrix structured into 6 parts
             *  _________________
             *  | * * * | * * * |
             *  | *(1)* | *(2)* |    _______
             *  |_*_*_*_|_*_*_*_|    | ,xx |
             *  | * * * | * * * |    | ,yy |
             *  | * * * | * * * |    | ,zz |
             *  | *(3)* | *(4)* |  * | ,yz |
             *  | * * * | * * * |    | ,xz |
             *  | * * * | * * * |    |_,xy_|
             *  |_*_*_*_|_*_*_*_|
             *  |_*(5)*_|_*(6)*_|
             */

            // Block (1) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3bt( 0, j ) = 3 * aJ2bt( 0, j ) * aJt( 0, j );
                aJ3bt( 1, j ) = 3 * aJ2bt( 1, j ) * aJt( 1, j );
                aJ3bt( 2, j ) = 3 * aJ2bt( 2, j ) * aJt( 2, j );
            }

            // Block (2) ------------------------------------------------
            aJ3bt( 0, 5 ) =
                    3 * aJ2bt( 0, 0 ) * aJt( 0, 1 ) +
                    3 * aJ2bt( 0, 1 ) * aJt( 0, 0 );
            aJ3bt( 1, 5 ) =
                    3 * aJ2bt( 1, 0 ) * aJt( 1, 1 ) +
                    3 * aJ2bt( 1, 1 ) * aJt( 1, 0 );
            aJ3bt( 2, 5 ) =
                    3 * aJ2bt( 2, 0 ) * aJt( 2, 1 ) +
                    3 * aJ2bt( 2, 1 ) * aJt( 2, 0 );

            aJ3bt( 0, 3 ) =
                    3 * aJ2bt( 0, 1 ) * aJt( 0, 2 ) +
                    3 * aJ2bt( 0, 2 ) * aJt( 0, 1 );
            aJ3bt( 1, 3 ) =
                    3 * aJ2bt( 1, 1 ) * aJt( 1, 2 ) +
                    3 * aJ2bt( 1, 2 ) * aJt( 1, 1 );
            aJ3bt( 2, 3 ) =
                    3 * aJ2bt( 2, 1 ) * aJt( 2, 2 ) +
                    3 * aJ2bt( 2, 2 ) * aJt( 2, 1 );

            aJ3bt( 0, 4 ) =
                    3 * aJ2bt( 0, 0 ) * aJt( 0, 2 ) +
                    3 * aJ2bt( 0, 2 ) * aJt( 0, 0 );
            aJ3bt( 1, 4 ) =
                    3 * aJ2bt( 1, 0 ) * aJt( 1, 2 ) +
                    3 * aJ2bt( 1, 2 ) * aJt( 1, 0 );
            aJ3bt( 2, 4 ) =
                    3 * aJ2bt( 2, 0 ) * aJt( 2, 2 ) +
                    3 * aJ2bt( 2, 2 ) * aJt( 2, 0 );

            // Block (3) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3bt( 3, j ) = 2 * aJ2bt( 5, j ) * aJt( 0, j )  +  aJ2bt( 0, j ) * aJt( 1, j );
                aJ3bt( 4, j ) = 2 * aJ2bt( 4, j ) * aJt( 0, j )  +  aJ2bt( 0, j ) * aJt( 2, j );
                aJ3bt( 5, j ) = 2 * aJ2bt( 5, j ) * aJt( 1, j )  +  aJ2bt( 1, j ) * aJt( 0, j );
                aJ3bt( 6, j ) = 2 * aJ2bt( 3, j ) * aJt( 1, j )  +  aJ2bt( 1, j ) * aJt( 2, j );
                aJ3bt( 7, j ) = 2 * aJ2bt( 4, j ) * aJt( 2, j )  +  aJ2bt( 2, j ) * aJt( 0, j );
                aJ3bt( 8, j ) = 2 * aJ2bt( 3, j ) * aJt( 2, j )  +  aJ2bt( 2, j ) * aJt( 1, j );
            }

            // Block (4) ------------------------------------------------
            aJ3bt( 3, 5 ) =  2 * aJ2bt( 5, 0 ) * aJt( 0, 1 )  +  2 * aJ2bt( 5, 1 ) * aJt( 0, 0 )  +  aJ2bt( 0, 1 ) * aJt( 1, 0 )  +  aJ2bt( 0, 0 ) * aJt( 1, 1 );
            aJ3bt( 4, 5 ) =  2 * aJ2bt( 4, 0 ) * aJt( 0, 1 )  +  2 * aJ2bt( 4, 1 ) * aJt( 0, 0 )  +  aJ2bt( 0, 1 ) * aJt( 2, 0 )  +  aJ2bt( 0, 0 ) * aJt( 2, 1 );
            aJ3bt( 5, 5 ) =  2 * aJ2bt( 5, 0 ) * aJt( 1, 1 )  +  2 * aJ2bt( 5, 1 ) * aJt( 1, 0 )  +  aJ2bt( 1, 1 ) * aJt( 0, 0 )  +  aJ2bt( 1, 0 ) * aJt( 0, 1 );
            aJ3bt( 6, 5 ) =  2 * aJ2bt( 3, 0 ) * aJt( 1, 1 )  +  2 * aJ2bt( 3, 1 ) * aJt( 1, 0 )  +  aJ2bt( 1, 1 ) * aJt( 2, 0 )  +  aJ2bt( 1, 0 ) * aJt( 2, 1 );
            aJ3bt( 7, 5 ) =  2 * aJ2bt( 4, 0 ) * aJt( 2, 1 )  +  2 * aJ2bt( 4, 1 ) * aJt( 2, 0 )  +  aJ2bt( 2, 1 ) * aJt( 0, 0 )  +  aJ2bt( 2, 0 ) * aJt( 0, 1 );
            aJ3bt( 8, 5 ) =  2 * aJ2bt( 3, 0 ) * aJt( 2, 1 )  +  2 * aJ2bt( 3, 1 ) * aJt( 2, 0 )  +  aJ2bt( 2, 1 ) * aJt( 1, 0 )  +  aJ2bt( 2, 0 ) * aJt( 1, 1 );

            aJ3bt( 3, 3 ) =  2 * aJ2bt( 5, 1 ) * aJt( 0, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 0, 1 )  +  aJ2bt( 0, 2 ) * aJt( 1, 1 )  +  aJ2bt( 0, 1 ) * aJt( 1, 2 );
            aJ3bt( 4, 3 ) =  2 * aJ2bt( 4, 1 ) * aJt( 0, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 0, 1 )  +  aJ2bt( 0, 2 ) * aJt( 2, 1 )  +  aJ2bt( 0, 1 ) * aJt( 2, 2 );
            aJ3bt( 5, 3 ) =  2 * aJ2bt( 5, 1 ) * aJt( 1, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 1, 1 )  +  aJ2bt( 1, 2 ) * aJt( 0, 1 )  +  aJ2bt( 1, 1 ) * aJt( 0, 2 );
            aJ3bt( 6, 3 ) =  2 * aJ2bt( 3, 1 ) * aJt( 1, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 1, 1 )  +  aJ2bt( 1, 2 ) * aJt( 2, 1 )  +  aJ2bt( 1, 1 ) * aJt( 2, 2 );
            aJ3bt( 7, 3 ) =  2 * aJ2bt( 4, 1 ) * aJt( 2, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 2, 1 )  +  aJ2bt( 2, 2 ) * aJt( 0, 1 )  +  aJ2bt( 2, 1 ) * aJt( 0, 2 );
            aJ3bt( 8, 3 ) =  2 * aJ2bt( 3, 1 ) * aJt( 2, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 2, 1 )  +  aJ2bt( 2, 2 ) * aJt( 1, 1 )  +  aJ2bt( 2, 1 ) * aJt( 1, 2 );

            aJ3bt( 3, 4 ) =  2 * aJ2bt( 5, 0 ) * aJt( 0, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 0, 0 )  +  aJ2bt( 0, 2 ) * aJt( 1, 0 )  +  aJ2bt( 0, 0 ) * aJt( 1, 2 );
            aJ3bt( 4, 4 ) =  2 * aJ2bt( 4, 0 ) * aJt( 0, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 0, 0 )  +  aJ2bt( 0, 2 ) * aJt( 2, 0 )  +  aJ2bt( 0, 0 ) * aJt( 2, 2 );
            aJ3bt( 5, 4 ) =  2 * aJ2bt( 5, 0 ) * aJt( 1, 2 )  +  2 * aJ2bt( 5, 2 ) * aJt( 1, 0 )  +  aJ2bt( 1, 2 ) * aJt( 0, 0 )  +  aJ2bt( 1, 0 ) * aJt( 0, 2 );
            aJ3bt( 6, 4 ) =  2 * aJ2bt( 3, 0 ) * aJt( 1, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 1, 0 )  +  aJ2bt( 1, 2 ) * aJt( 2, 0 )  +  aJ2bt( 1, 0 ) * aJt( 2, 2 );
            aJ3bt( 7, 4 ) =  2 * aJ2bt( 4, 0 ) * aJt( 2, 2 )  +  2 * aJ2bt( 4, 2 ) * aJt( 2, 0 )  +  aJ2bt( 2, 2 ) * aJt( 0, 0 )  +  aJ2bt( 2, 0 ) * aJt( 0, 2 );
            aJ3bt( 8, 4 ) =  2 * aJ2bt( 3, 0 ) * aJt( 2, 2 )  +  2 * aJ2bt( 3, 2 ) * aJt( 2, 0 )  +  aJ2bt( 2, 2 ) * aJt( 1, 0 )  +  aJ2bt( 2, 0 ) * aJt( 1, 2 );

            // Block (5) ------------------------------------------------
            for( uint j=0; j<3; ++j )
            {
                aJ3bt( 9, j ) =
                        aJ2bt( 4, j ) * aJt( 1, j ) +
                        aJ2bt( 3, j ) * aJt( 0, j ) +
                        aJ2bt( 5, j ) * aJt( 2, j );
            }

            // Block (6) ------------------------------------------------
            aJ3bt( 9, 5 ) =
                    aJ2bt( 5, 0 ) * aJt( 2, 1 ) + aJ2bt( 5, 1 ) * aJt( 2, 0 ) +
                    aJ2bt( 3, 0 ) * aJt( 0, 1 ) + aJ2bt( 3, 1 ) * aJt( 0, 0 ) +
                    aJ2bt( 4, 0 ) * aJt( 1, 1 ) + aJ2bt( 4, 1 ) * aJt( 1, 0 );

            aJ3bt( 9, 3 ) =
                    aJ2bt( 5, 1 ) * aJt( 2, 2 ) + aJ2bt( 5, 2 ) * aJt( 2, 1 ) +
                    aJ2bt( 3, 1 ) * aJt( 0, 2 ) + aJ2bt( 3, 2 ) * aJt( 0, 1 ) +
                    aJ2bt( 4, 1 ) * aJt( 1, 2 ) + aJ2bt( 4, 2 ) * aJt( 1, 1 );

            aJ3bt( 9, 4 ) =
                    aJ2bt( 5, 0 ) * aJt( 2, 2 ) + aJ2bt( 5, 2 ) * aJt( 2, 0 ) +
                    aJ2bt( 3, 0 ) * aJt( 0, 2 ) + aJ2bt( 3, 2 ) * aJt( 0, 0 ) +
                    aJ2bt( 4, 0 ) * aJt( 1, 2 ) + aJ2bt( 4, 2 ) * aJt( 1, 0 );

            // third help matrix
            aJ3ct = ad3NdXi3 * aXHat;
        }

        //------------------------------------------------------------------------------

        void Geometry_Interpolator::set_function_pointers()
        {
            // get number of dimensions and set pointer to function
            // for second space derivative
            switch ( mNumSpaceDim )
            {
                case 1 :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_1d;
                    mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_1d;
                    break;
                }
                case 2 :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_2d;
                    mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_2d;
                    break;
                }
                case 3 :
                {
                    mSecondDerivativeMatricesSpace = this->eval_matrices_for_second_derivative_3d;
                    mThirdDerivativeMatricesSpace  = this->eval_matrices_for_third_derivative_3d;
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
                case 1 :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_1d;
                    break;
                }
                case 2 :
                {
                    mSecondDerivativeMatricesTime = this->eval_matrices_for_second_derivative_2d;
                    break;
                }
                case 3 :
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

            // switch on geometry type
            if( mSpaceSideset )
            {
                switch( mGeometryType )
                {
                    case mtk::Geometry_Type::LINE :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_side_line;
                        mNormalFunc    = &Geometry_Interpolator::eval_normal_side_line;
                        break;
                    }
                    case mtk::Geometry_Type::TRI :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_side_tri;
                        mNormalFunc    = &Geometry_Interpolator::eval_normal_side_tri;
                        break;
                    }
                    case mtk::Geometry_Type::QUAD :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_side_quad;
                        mNormalFunc    = &Geometry_Interpolator::eval_normal_side_quad;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented side space geometry type ");
                        break;
                    }
                }
            }
            else
            {
                switch( mGeometryType )
                {
                    case mtk::Geometry_Type::LINE :
                    case mtk::Geometry_Type::QUAD :
                    case mtk::Geometry_Type::HEX :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_line_quad_hex;
                        break;
                    }

                    case mtk::Geometry_Type::TRI :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_tri;
                        break;
                    }

                    case mtk::Geometry_Type::TET :
                    {
                        mSpaceDetJFunc = &Geometry_Interpolator::eval_space_detJ_bulk_tet;
                        break;
                    }

                    default :
                    {
                        MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented space geometry type ");
                        break;
                    }
                }
            }

            // switch Geometry_Type
            switch ( mTimeGeometryType )
            {
                case mtk::Geometry_Type::POINT :
                {
                    mTimeDetJFunc = &Geometry_Interpolator::eval_time_detJ_side;
                    break;
                }
                case mtk::Geometry_Type::LINE :
                {
                    mTimeDetJFunc = &Geometry_Interpolator::eval_time_detJ_bulk;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Geometry_Interpolator::set_function_pointers - unknown or not implemented time geometry type ");
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------

        moris::Cell< moris::Cell < moris_index > > Geometry_Interpolator::get_face_vertices_ordinals(
                uint aNumSpaceBases )
        {
            // init the vertices per facet
            moris::Cell< moris::Cell < moris_index > > tVerticesPerFace;

            // depending on the parent geometry
            switch ( mGeometryType )
            {
                case mtk::Geometry_Type::LINE :
                {
                    switch ( aNumSpaceBases )
                    {
                        case 1 :
                            tVerticesPerFace = { { 0 } };
                            break;
                        case 2 :
                            tVerticesPerFace = { { 0 }, { 1 } };
                            break;
                        case 3 :
                            tVerticesPerFace = { { 0 }, { 1 }, { 2 } };
                            break;
                        case 4 :
                            tVerticesPerFace = { { 0 }, { 1 }, { 2 }, { 3 } };
                            break;
                        default:
                            MORIS_ERROR( false, "Geometry_Interpolator::get_face_vertices_ordinals - LINE order not implemented " );
                            break;
                    }
                    break;
                }

                case mtk::Geometry_Type::QUAD :
                {
                    switch ( aNumSpaceBases )
                    {
                        case 4 :
                            tVerticesPerFace = {
                                    { 0, 1 },
                                    { 1, 2 },
                                    { 2, 3 },
                                    { 3, 0 } };
                            break;
                        case 8 :
                            tVerticesPerFace = {
                                    { 0, 1, 4 },
                                    { 1, 2, 5 },
                                    { 2, 3, 6 },
                                    { 3, 0, 7 } };
                            break;
                        case 9 :
                            tVerticesPerFace = {
                                    { 0, 1, 4 },
                                    { 1, 2, 5 },
                                    { 2, 3, 6 },
                                    { 3, 0, 7 } };
                            break;
                        case 16 :
                            tVerticesPerFace = {
                                    { 0, 1,  4,  5 },
                                    { 1, 2,  6,  7 },
                                    { 2, 3,  8,  9 },
                                    { 3, 0, 10, 11 } };
                            break;
                        default:
                            MORIS_ERROR( false, "Geometry_Interpolator::get_face_vertices_ordinals - QUAD order not implemented " );
                            break;
                    }
                    break;
                }

                case mtk::Geometry_Type::HEX :
                {
                    switch( aNumSpaceBases )
                    {
                        case 8 :
                            tVerticesPerFace = {
                                    { 0, 1, 5, 4 },
                                    { 1, 2, 6, 5 },
                                    { 2, 3, 7, 6 },
                                    { 0, 4, 7, 3 },
                                    { 0, 3, 2, 1 },
                                    { 4, 5, 6, 7 } };
                            break;
                        case 20 :
                            tVerticesPerFace = {
                                    { 0, 1, 5, 4,  8, 13, 16, 12 },
                                    { 1, 2, 6, 5,  9, 14, 17, 13 },
                                    { 2, 3, 7, 6, 10, 15, 18, 14 },
                                    { 0, 4, 7, 3, 12, 19, 15, 11 },
                                    { 0, 3, 2, 1, 11, 10,  9,  8 },
                                    { 4, 5, 6, 7, 16, 17, 18, 19 } };
                            break;
                        case 27 :
                            tVerticesPerFace = {
                                    { 0, 1, 5, 4,  8, 13, 16, 12, 23 },
                                    { 1, 2, 6, 5,  9, 14, 17, 13, 24 },
                                    { 2, 3, 7, 6, 10, 15, 18, 14, 26 },
                                    { 0, 4, 7, 3, 12, 19, 15, 11, 23 },
                                    { 0, 3, 2, 1, 11, 10,  9,  8, 21 },
                                    { 4, 5, 6, 7, 16, 17, 18, 19, 22 } };
                            break;
                        case 64 :
                            tVerticesPerFace = {
                                    { 0, 1, 5, 4,  8,  9, 16, 17, 25, 24, 13, 12, 36, 37, 38, 39 },
                                    { 1, 2, 6, 5, 14, 15, 20, 21, 29, 28, 17, 16, 44, 45, 46, 47 },
                                    { 2, 3, 7, 6, 18, 19, 22, 23, 31, 30, 21, 20, 48, 49, 50, 51 },
                                    { 0, 4, 7, 3, 12, 13, 26, 27, 23, 22, 11, 10, 40, 41, 42, 43 },
                                    { 0, 3, 2, 1, 10, 11, 19, 18, 15, 14,  9,  8, 32, 33, 34, 35 },
                                    { 4, 5, 6, 7, 24, 25, 28, 29, 30, 31, 27, 26, 52, 53, 54, 55 }};
                            break;
                        default:
                            MORIS_ERROR( false, "Geometry_Interpolator::get_face_vertices_ordinals - HEX order not implemented " );
                            break;
                    }
                    break;
                }

                case mtk::Geometry_Type::TRI :
                {
                    switch( aNumSpaceBases )
                    {
                        case 3 :
                            tVerticesPerFace = {
                                    { 0, 1 },
                                    { 1, 2 },
                                    { 2, 0 }};
                            break;
                        case 6 :
                            tVerticesPerFace = {
                                    { 0, 1, 3 },
                                    { 1, 2, 4 },
                                    { 2, 0, 5 }};
                            break;
                        case 10 :
                            tVerticesPerFace = {
                                    { 0, 1, 3, 4 },
                                    { 1, 2, 5, 6 },
                                    { 2, 0, 7, 8 }};
                            break;
                        default:
                            MORIS_ERROR( false, "Geometry_Interpolator::get_face_vertices_ordinals - TRI order not implemented " );
                            break;
                    }
                    break;
                }

                case mtk::Geometry_Type::TET :
                {
                    switch( aNumSpaceBases )
                    {
                        case 4 :
                            tVerticesPerFace = {
                                    { 0, 1, 3 },
                                    { 1, 2, 3 },
                                    { 0, 3, 2 },
                                    { 0, 2, 1 }};
                            break;
                        case 10 :
                            tVerticesPerFace = {
                                    { 0, 1, 3, 4, 8, 7 },
                                    { 1, 2, 3, 5, 9, 8 },
                                    { 0, 3, 2, 7, 9, 6 },
                                    { 0, 2, 1, 6, 5, 4 }};
                            break;
                        case 20 :
                            tVerticesPerFace = {
                                    { 0, 1, 3,  4,  5, 12, 13, 11, 10, 17 },
                                    { 1, 2, 3,  6,  7, 14, 15, 13, 12, 18 },
                                    { 0, 3, 2, 10, 11, 15, 14,  9,  8, 19 },
                                    { 0, 2, 1,  8,  9,  7,  6,  5,  4, 16 }};
                            break;
                        default:
                            MORIS_ERROR( false, "Geometry_Interpolator::get_face_vertices_ordinals - TET order not implemented " );
                            break;
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "Geometry_Interpolator::get_space_sideset_param_coeff - undefined geometry type " );
                    break;
                }
            }
            return tVerticesPerFace;
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
