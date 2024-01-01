/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_FEM_IWG_Compressible_NS_Bulk_Stabilization.cpp  
 * 
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG_Compressible_NS_Bulk.hpp"
#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_inv.hpp"
#include "fn_sqrtmat.hpp"
#include "fn_sylvester.hpp"

// debug - output to hdf5
#include "fn_max.hpp"
#include "paths.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::A0inv()
        {
            // check if the operator has already been evaluated
            if ( !mA0invEval )
            {
                return mA0inv;
            }

            // update eval flag
            mA0invEval = false;

            // get A0 for the given variable set and invert it
            mA0inv = inv( this->A( 0 ) );

            // return
            return mA0inv;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::G()
        {
            // check if G has already been evaluated
            if ( !mGEval )
            {
                return mG;
            }

            // update eval flag
            mGEval = false;
            
            // get the space jacobian from IP geometry interpolator
            const Matrix< DDRMat > & tInvSpaceJacobian =
                    mLeaderFIManager->get_IP_geometry_interpolator()->inverse_space_jacobian();

            // compute G based on number of spatial dimensions
            if ( this->num_space_dims() == 2 )
            {
                // set size for mG
                mG.set_size( 2, 2);

                // fill aGij = sum_d dxi_d/dx_i dxi_d/dx_j
                mG( 0, 0 ) = std::pow( tInvSpaceJacobian( 0, 0 ), 2.0 ) + std::pow( tInvSpaceJacobian( 0, 1 ), 2.0 );
                mG( 0, 1 ) = tInvSpaceJacobian( 0, 0 ) * tInvSpaceJacobian( 1, 0 ) + tInvSpaceJacobian( 0, 1 ) * tInvSpaceJacobian( 1, 1 );
                mG( 1, 0 ) = mG( 0, 1 );
                mG( 1, 1 ) = std::pow( tInvSpaceJacobian( 1, 0 ), 2.0 ) + std::pow( tInvSpaceJacobian( 1, 1 ), 2.0 );
            }
            else if ( this->num_space_dims() == 3 ) 
            {
                // set size for mG
                mG.set_size( 3, 3);

                // fill aGij = sum_d dxi_d/dx_i dxi_d/dx_j
                mG( 0, 0 ) = std::pow( tInvSpaceJacobian( 0, 0 ), 2.0 )
                + std::pow( tInvSpaceJacobian( 0, 1 ), 2.0 )
                + std::pow( tInvSpaceJacobian( 0, 2 ), 2.0 );
                mG( 0, 1 ) = tInvSpaceJacobian( 0, 0 ) * tInvSpaceJacobian( 1, 0 )
                                + tInvSpaceJacobian( 0, 1 ) * tInvSpaceJacobian( 1, 1 )
                                + tInvSpaceJacobian( 0, 2 ) * tInvSpaceJacobian( 1, 2 );
                mG( 0, 2 ) = tInvSpaceJacobian( 0, 0 ) * tInvSpaceJacobian( 2, 0 )
                                + tInvSpaceJacobian( 0, 1 ) * tInvSpaceJacobian( 2, 1 )
                                + tInvSpaceJacobian( 0, 2 ) * tInvSpaceJacobian( 2, 2 );
                mG( 1, 0 ) = mG( 0, 1 );
                mG( 1, 1 ) = std::pow( tInvSpaceJacobian( 1, 0 ), 2.0 )
                                + std::pow( tInvSpaceJacobian( 1, 1 ), 2.0 )
                                + std::pow( tInvSpaceJacobian( 1, 2 ), 2.0 );
                mG( 1, 2 ) = tInvSpaceJacobian( 1, 0 ) * tInvSpaceJacobian( 2, 0 )
                                + tInvSpaceJacobian( 1, 1 ) * tInvSpaceJacobian( 2, 1 )
                                + tInvSpaceJacobian( 1, 2 ) * tInvSpaceJacobian( 2, 2 );
                mG( 2, 0 ) = mG( 0, 2 );
                mG( 2, 1 ) = mG( 1, 2 );
                mG( 2, 2 ) = std::pow( tInvSpaceJacobian( 2, 0 ), 2.0 )
                                + std::pow( tInvSpaceJacobian( 2, 1 ), 2.0 )
                                + std::pow( tInvSpaceJacobian( 2, 2 ), 2.0 );
            } 
            else
            {
                MORIS_ERROR( false, "IWG_Compressible_NS_Bulk::G() - Number of spatial dimensions must be 2 or 3." );
            }

            // return mG
            return mG;  
        }  

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::M()
        {
            // check if the operator has already been evaluated
            if ( !mMEval )
            {
                return mM;
            }

            // update eval flag
            mMEval = false;

            // compute M ....

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // FIXME: C-operator for body forces and heat load ignored for right now

            // get the time jacobian from IP geometry interpolator
            real tdTaudt =
                    mLeaderFIManager->get_IP_geometry_interpolator()->inverse_time_jacobian()( 0 );

            // get identity matrix of correct size
            Matrix< DDRMat > tIdentity;
            eye( tNumStateVars, tNumStateVars, tIdentity );

            // initialize M with time scaling term
            mM = tdTaudt * tdTaudt * tIdentity;

            // get subview of mM for += operatorions
            auto tM = mM( { 0, tNumStateVars - 1 }, { 0, tNumStateVars - 1 } );

            // add body force term to M
            Matrix< DDRMat > tC = this->C() * this->A0inv();
            tM += tC * tC;

            // add loop over A and K terms
            for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
            {
                for ( uint kDim = 0; kDim < this->num_space_dims(); kDim++ )
                {
                    // add contribution from the A-terms
                    tM += this->G()( jDim, kDim ) * this->A( jDim + 1 ) * this->A0inv() * this->A( kDim + 1 ) * this->A0inv();

                    for ( uint lDim = 0; lDim < this->num_space_dims(); lDim++ )
                    {
                        for ( uint mDim = 0; mDim < this->num_space_dims(); mDim++ )
                        {
                            // clang-format off
                            // add contribution from the K-terms
                            tM += this->G()( kDim, jDim ) * this->G()( lDim, mDim ) * 
                                    this->K( kDim, lDim ) * this->A0inv() *
                                    this->K( mDim, jDim ) * this->A0inv(); 
                            // clang-format on
                        }
                    }
                }
            }

            // return value
            return mM;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::Minv()
        {
            // check if the operator has already been evaluated
            if ( !mMinvEval )
            {
                return mMinv;
            }

            // update eval flag
            mMinvEval = false;

            // get A0 for the given variable set and invert it
            mMinv = inv( this->M() );

            // return
            return mMinv;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::SqrtMinv()
        {
            // check if the operator has already been evaluated
            if ( !mSqrtMinvEval )
            {
                return mSqrtMinv;
            }

            // update eval flag
            mSqrtMinvEval = false;

            // get the squareroot of the inverse of M
            moris::sqrtmat( this->Minv(), mSqrtMinv );

            // return
            return mSqrtMinv;
        }        

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::Tau()
        {
            // check if Tau has already been evaluated
            if ( !mTauEval )
            {
                return mTau;
            }

            // update eval flag
            mTauEval = false;

            // compute Tau
            mTau = this->A0inv() * this->SqrtMinv();

            // return
            return mTau;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::dMdY( const uint aYind )
        {
            // check if Tau has already been evaluated
            if ( !mdMdYEval )
            {
                return mdMdY( aYind );
            }

            // update eval flag
            mdMdYEval = false;

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // initialize cell for storage
            Matrix< DDRMat > tZeroMatrix( tNumStateVars, tNumStateVars, 0.0 );
            mdMdY.assign( tNumStateVars, tZeroMatrix );

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropMu = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // for each state variable compute the derivative
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // get subview for += operations
                auto tdMdVar = mdMdY( iVar )( { 0, tNumStateVars - 1 }, { 0, tNumStateVars - 1 } );

                // add body force term to M
                Matrix< DDRMat > tdCdY = this->dCdY(iVar) * this->A0inv() + this->C() * this->dA0invdY( iVar );
                Matrix< DDRMat > tC    = this->C() * this->A0inv();
                tdMdVar += tdCdY * tC + tC * tdCdY;

                // get the variable derivs for the K matrices
                moris::Vector< moris::Vector< Matrix< DDRMat > > > tdKdY;
                eval_dKdY( tPropMu, tPropKappa, mLeaderFIManager, iVar, tdKdY );

                // loops for addition over indices j,k,l,m
                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    // get the variable derivs for the A matrices
                    Matrix< DDRMat > tdAjdY;
                    eval_dAdY( tMM, tCM, mLeaderFIManager, mResidualDofType, jDim + 1, iVar, tdAjdY );

                    for ( uint kDim = 0; kDim < this->num_space_dims(); kDim++ )
                    {
                        // get the variable derivs for the A matrices
                        Matrix< DDRMat > tdAkdY;
                        eval_dAdY( tMM, tCM, mLeaderFIManager, mResidualDofType, kDim + 1, iVar, tdAkdY );

                        // clang-format off
                        // add contribution from the A-terms
                        tdMdVar += this->G()( jDim, kDim ) * (  // tdMdVar
                                             tdAjdY * this->A0inv()          * this->A( kDim + 1 ) * this->A0inv() + 
                                this->A( jDim + 1 ) * this->dA0invdY( iVar ) * this->A( kDim + 1 ) * this->A0inv() + 
                                this->A( jDim + 1 ) * this->A0inv()          *              tdAkdY * this->A0inv() +
                                this->A( jDim + 1 ) * this->A0inv()          * this->A( kDim + 1 ) * this->dA0invdY( iVar ) );

                        for ( uint lDim = 0; lDim < this->num_space_dims(); lDim++ )
                        {
                            for ( uint mDim = 0; mDim < this->num_space_dims(); mDim++ )
                            {
                                // add contribution from the K-terms
                                tdMdVar += this->G()( kDim, jDim ) * this->G()( lDim, mDim ) * ( // tdMdVar
                                        tdKdY( kDim )( lDim ) * this->A0inv()          * this->K( mDim, jDim ) * this->A0inv() + 
                                        this->K( kDim, lDim ) * this->dA0invdY( iVar ) * this->K( mDim, jDim ) * this->A0inv() + 
                                        this->K( kDim, lDim ) * this->A0inv()          * tdKdY( mDim )( jDim ) * this->A0inv() +
                                        this->K( kDim, lDim ) * this->A0inv()          * this->K( mDim, jDim ) * this->dA0invdY( iVar ) ); 
                            }
                        }
                        // clang-format on
                    }
                }
            }

            // return
            return mdMdY( aYind );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::dA0invdY( const uint aYind )
        {
            // check if Tau has already been evaluated
            if ( !mdA0invdYEval )
            {
                return mdA0invdY( aYind );
            }

            // update eval flag
            mdA0invdYEval = false;

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;

            // initialize
            mdA0invdY.resize( tNumStateVars );

            // for each state variable compute the derivative
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // get the variable derivs for the A matrices
                Matrix< DDRMat > tdA0dY;
                eval_dAdY( tMM, tCM, mLeaderFIManager, mResidualDofType, 0, iVar, tdA0dY );

                // compute the state var deriv of the inverse of A0
                mdA0invdY( iVar ) = -1.0 * this->A0inv() * tdA0dY * this->A0inv();
            }

            // return
            return mdA0invdY( aYind );
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::dTaudY( const Matrix< DDRMat > aVR )
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get number of state variables
            uint tNumStateVars = this->num_space_dims() + 2;
            
            // -------------------------------------------------
            // STEP 1: get the state variable derivatives of M
            // see two functions above

            // -------------------------------------------------
            // STEP 2: get the state variable derivatives of the inverse of M

            // initialize
            moris::Vector< Matrix< DDRMat > > tMinusdMinvdY( tNumStateVars );
            
            // for each state variable compute the derivative
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // compute the state var deriv
                tMinusdMinvdY( iVar ) = this->Minv() * this->dMdY( iVar ) * this->Minv();
            }

            // -------------------------------------------------
            // STEP 3: get the square root of the state variable derivatives 
            //         of the inverse of M using the Sylvester equation

            // initialize
            moris::Vector< Matrix< DDRMat > > tdSqrtMinvdY( tNumStateVars );
            
            // for each state variable compute the derivative
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // solve Sylvester equation and store solutions for each variable derivative
                sylvester( this->SqrtMinv(), this->SqrtMinv(), tMinusdMinvdY( iVar ), tdSqrtMinvdY( iVar ) );
            }

            // -------------------------------------------------
            // STEP 4: get the state variable derivatives of the inverse of A0
            // see function above

            // -------------------------------------------------
            // STEP 5: post-multiplication with the input vector

            // initialize
            mdTaudY.set_size( tNumStateVars, tNumStateVars, 0.0 );

            // for each state variable compute the derivative
            for ( uint iVar = 0; iVar < tNumStateVars; iVar++ )
            {
                // perform multiplication and put everything in
                mdTaudY( { 0, tNumStateVars - 1 }, { iVar, iVar } ) = 
                        this->A0inv() * tdSqrtMinvdY( iVar ) * aVR +
                        this->dA0invdY( iVar ) * this->SqrtMinv() * aVR;
            }

            // -------------------------------------------------
            // return value
            return mdTaudY;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::LY()
        {
            // check if LY already been evaluated
            if ( !mLYEval )
            {
                return mLY;
            }

            // set the eval flag
            mLYEval = false;  

            // initialize LY with A0 term
            mLY = this->A( 0 ) * this->dYdt();

            // get subview for += operations
            auto tLY = mLY( { 0, mLY.n_rows() - 1 }, { 0, mLY.n_cols() - 1 } );

            // add contribution due to body loads
            tLY += this->C() * this->Y();

            // loop over A and K matrices
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                tLY += ( this->A( iDim + 1 ) - this->Kiji( iDim ) ) * this->dYdx( iDim ); 

                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    tLY -= this->K( iDim, jDim ) * this->d2Ydx2( iDim, jDim );
                }
            }

            // return value
            return mLY;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::dLdDofY()
        {
            // check if LY already been evaluated
            if ( !mLDofYEval )
            {
                return mdLdDofY;
            }

            // set the eval flag
            mLDofYEval = false;  

            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropMu = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // initialize cell containing A-matrices pre-multiplied with the state variable vector
            moris::Vector< Matrix< DDRMat > > tdAjdY_Yj( this->num_space_dims() + 1 );

            // initialize cell containing Kij,i-matrices pre-multiplied with the state variable vector
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tdKijidY_Yj( this->num_space_dims() );

            // initialize cell containing Kij,i-matrices pre-multiplied with the state variable vector
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tdKijdY_Yij( this->num_space_dims() );
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++)
            {
                tdKijdY_Yij( iDim ).resize( this->num_space_dims() );
            }

            // get dA0/dY * Y,t
            eval_dAdY_VR( tMM, tCM, mLeaderFIManager, mResidualDofType, this->dYdt(), 0, tdAjdY_Yj( 0 ) );

            // compute A(0) term
            mdLdDofY = tdAjdY_Yj( 0 ) * this->W();

            // get subview of matrix for += operations
            auto tdLdDofY = mdLdDofY( { 0, mdLdDofY.n_rows() - 1 }, { 0, mdLdDofY.n_cols() - 1 } );

            // add contribution due to body loads
            tdLdDofY += this->dCdY_VR( this->Y() ) * this->W();

            // go over all Aj*Y,j and Kij,i*Y,j and Kij*Y,ij terms and add up
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                // get dAj/dY * Y,j
                eval_dAdY_VR( tMM, tCM, mLeaderFIManager, mResidualDofType, this->dYdx( iDim ), iDim + 1, tdAjdY_Yj( iDim + 1 ) );

                // add contributions from A-matrices
                tdLdDofY += tdAjdY_Yj( iDim + 1 )  * this->W();

                // get dKij,i/dY * Y,j
                eval_dKijidY_VR( tPropMu, tPropKappa, mLeaderFIManager, this->dYdx( iDim ), iDim, tdKijidY_Yj( iDim ) );

                // add contributions from Kij,i-matrices
                // tdLdDofY -= tdKijidY_Yj( iDim )( 0 ) * this->W();

                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    // add contributions from Kij,i-matrices
                    tdLdDofY -= tdKijidY_Yj( iDim )( jDim + 1 ) * this->dWdx( jDim );

                    // get dKij/dY * Y,ij
                    eval_dKdY_VR( tPropMu, tPropKappa, mLeaderFIManager, this->d2Ydx2( iDim, jDim ), iDim, jDim, tdKijdY_Yij( iDim )( jDim ) );
                    
                    // add contributions from K-matrices
                    tdLdDofY -= tdKijdY_Yij( iDim )( jDim ) * this->W();
                }
            }

            // return value
            return mdLdDofY;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::LW()
        {
            // check if LY already been evaluated
            if ( !mLWEval )
            {
                return mLW;
            }

            // set the eval flag
            mLWEval = false;  

            // initialize LW with A0 term
            mLW = this->A( 0 ) * this->dWdt();

            // get subview for += operations
            auto tLW = mLW( { 0, mLW.n_rows() - 1 }, { 0, mLW.n_cols() - 1 } );

            // add contribution due to body loads
            tLW += this->C() * this->W();

            // loop over A- and K-matrices
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                tLW += ( this->A( iDim + 1 ) - this->Kiji( iDim ) ) * this->dWdx( iDim ); 

                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    tLW -= this->K( iDim, jDim ) * this->d2Wdx2( iDim, jDim );
                }
            }

            // return value
            return mLW;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::GLSTestFunc()
        {
            // check if LY already been evaluated
            if ( !mGLSTestFuncEval )
            {
                return mGLSTestFunc;
            }

            // set the eval flag
            mGLSTestFuncEval = false;  

            // initialize LW with A0 term
            mGLSTestFunc = this->dWdt_trans() * this->A( 0 );

            // get subview for += operations
            auto tGLSTF = mGLSTestFunc( { 0, mGLSTestFunc.n_rows() - 1 }, { 0, mGLSTestFunc.n_cols() - 1 } );

            // loop over A- and K-matrices
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                tGLSTF += this->dWdx_trans( iDim ) * ( this->A( iDim + 1 ) - this->Kiji( iDim ) ); 

                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    tGLSTF -= trans( this->d2Wdx2( iDim, jDim ) ) * this->K( iDim, jDim );
                }
            }

            // return value
            return mGLSTestFunc;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::dLdDofW( const Matrix< DDRMat > & aVL )
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropMu = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // initialize cell containing A-matrices pre-multiplied with VL
            moris::Vector< Matrix< DDRMat > > tVLdAjdY( this->num_space_dims() + 1 );

            // initialize cell containing Kij,i-matrices pre-multiplied with VL
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tVLdKijidY( this->num_space_dims() );

            // initialize cell containing Kij,i-matrices pre-multiplied with VL
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tVLdKijdY( this->num_space_dims() );
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++)
            {
                tVLdKijdY( iDim ).resize( this->num_space_dims() );
            }

            // get VL * dA0/dY
            eval_VL_dAdY( tMM, tCM, mLeaderFIManager, mResidualDofType, aVL, 0, tVLdAjdY( 0 ) );

            // compute A(0) term
            mdLdDofW = this->dWdt_trans() * tVLdAjdY( 0 ) * this->W();

            // get subview of matrix for += operations
            auto tdLdDofW = mdLdDofW( { 0, mdLdDofW.n_rows() - 1 }, { 0, mdLdDofW.n_cols() - 1 } );

            // add contribution due to body loads // FIXME: dCdY_VR can only be applied because C is symmetric, dCdY_VL is needed
            tdLdDofW += this->W_trans() * this->dCdY_VR( aVL ) * this->W();

            // go over all VL*Aj and VL*Kij,i and VL*Kij terms and add up
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                // get VL * dAj/dY
                eval_VL_dAdY( tMM, tCM, mLeaderFIManager, mResidualDofType, aVL, iDim + 1, tVLdAjdY( iDim + 1 ) );

                // add contributions from A-matrices
                tdLdDofW += this->dWdx_trans( iDim ) * tVLdAjdY( iDim + 1 )  * this->W();

                // get VL * dKij,i/dY
                eval_VL_dKijidY( tPropMu, tPropKappa, mLeaderFIManager, aVL, iDim, tVLdKijidY( iDim ) );

                // add contributions from Kij,i-matrices
                // tdLdDofW -= this->dWdx_trans( iDim ) * tVLdKijidY( iDim )( 0 ) * this->W();

                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    // add contributions from Kij,i-matrices
                    tdLdDofW -= this->dWdx_trans( iDim ) * tVLdKijidY( iDim )( jDim + 1 ) * this->dWdx( jDim );

                    // get VL * dKij/dY
                    eval_VL_dKdY( tPropMu, tPropKappa, mLeaderFIManager, aVL, iDim, jDim, tVLdKijdY( iDim )( jDim ) );

                    // add contributions from K-matrices
                    tdLdDofW -= trans( this->d2Wdx2( iDim, jDim ) ) * tVLdKijdY( iDim )( jDim ) * this->W();
                }
            }

            // return value
            return mdLdDofW;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & IWG_Compressible_NS_Bulk::dGLSTestFuncdDof( const Matrix< DDRMat > & aVR )
        {
            // get the material and constitutive models
            std::shared_ptr< Material_Model > tMM = mLeaderMM( static_cast< uint >( IWG_Material_Type::FLUID_MM ) );
            std::shared_ptr< Constitutive_Model > tCM = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_CM ) );

            // get the properties
            std::shared_ptr< Property > tPropMu = mLeaderProp( static_cast< uint >( IWG_Property_Type::DYNAMIC_VISCOSITY ) );
            std::shared_ptr< Property > tPropKappa = mLeaderProp( static_cast< uint >( IWG_Property_Type::THERMAL_CONDUCTIVITY ) );

            // initialize cell containing A-matrices pre-multiplied with VR
            moris::Vector< Matrix< DDRMat > > tdAjdYVR( this->num_space_dims() + 1 );

            // initialize cell containing Kij,i-matrices pre-multiplied with VR
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tdKijidYVR( this->num_space_dims() );

            // initialize cell containing Kij,i-matrices pre-multiplied with VL
            moris::Vector< moris::Vector< Matrix< DDRMat > > > tdKijdYVR( this->num_space_dims() );
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++)
            {
                tdKijdYVR( iDim ).resize( this->num_space_dims() );
            }

            // get dA0/dY * VR
            eval_dAdY_VR( tMM, tCM, mLeaderFIManager, mResidualDofType, aVR, 0, tdAjdYVR( 0 ) );

            // compute A(0) term
            mdGLSTestFuncdDof = this->dWdt_trans() * tdAjdYVR( 0 ) * this->W();

            // get subview of matrix for += operations
            auto tdGLSTFdDof = mdGLSTestFuncdDof( { 0, mdGLSTestFuncdDof.n_rows() - 1 }, { 0, mdGLSTestFuncdDof.n_cols() - 1 } );

            // go over all Aj*VR and Kij,i*VR and Kij*VR terms and add up
            for ( uint iDim = 0; iDim < this->num_space_dims(); iDim++ )
            {
                // get dAj/dY * VR
                eval_dAdY_VR( tMM, tCM, mLeaderFIManager, mResidualDofType, aVR, iDim + 1, tdAjdYVR( iDim + 1 ) );

                // add contributions from A-matrices
                tdGLSTFdDof += this->dWdx_trans( iDim ) * tdAjdYVR( iDim + 1 )  * this->W();

                // get dKij,i/dY * VR
                eval_dKijidY_VR( tPropMu, tPropKappa, mLeaderFIManager, aVR, iDim, tdKijidYVR( iDim ) );

                // add contributions from Kij,i-matrices
                // tdGLSTFdDof -= this->dWdx_trans( iDim ) * tdKijidYVR( iDim )( 0 ) * this->W();

                for ( uint jDim = 0; jDim < this->num_space_dims(); jDim++ )
                {
                    // add contributions from Kij,i-matrices
                    tdGLSTFdDof -= this->dWdx_trans( iDim ) * tdKijidYVR( iDim )( jDim + 1 ) * this->dWdx( jDim );

                    // get dKij/dY * VR
                    eval_dKdY_VR( tPropMu, tPropKappa, mLeaderFIManager, aVR, iDim, jDim, tdKijdYVR( iDim )( jDim ) );

                    // add contributions from K-matrices
                    tdGLSTFdDof -= trans( this->d2Wdx2( iDim, jDim ) ) * tdKijdYVR( iDim )( jDim ) * this->W();
                }
            }

            // return value
            return mdGLSTestFuncdDof;
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
