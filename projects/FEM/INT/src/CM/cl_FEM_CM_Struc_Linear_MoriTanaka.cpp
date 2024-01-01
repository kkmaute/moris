/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear_MoriTanaka.cpp
 *
 */

#include "cl_FEM_CM_Struc_Linear_MoriTanaka.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        CM_Struc_Linear_MoriTanaka::CM_Struc_Linear_MoriTanaka()
                : mEshelbyTensor( 6, 6, 0.0 )
                , mConstMatrix( 6, 6, 0.0 )
                , mConstFiber( 6, 6, 0.0 )
        {
            //  set the property pointer cell size
            mProperties.resize( mProperties.size() + static_cast< uint >( CM_Property_Type_MT::MAX_ENUM ), nullptr );

            // This index is due to the existence of properties in the parent class
            uint tCurrentIndexOffSet = static_cast< uint >( CM_Property_Type_Lin::MAX_ENUM );

            // populate the map
            mPropertyMap[ "YoungsModulusMatrix" ] = static_cast< uint >( CM_Property_Type_MT::EMOD1 ) + tCurrentIndexOffSet;
            mPropertyMap[ "PoissonRatioMatrix" ]  = static_cast< uint >( CM_Property_Type_MT::NU1 ) + tCurrentIndexOffSet;
            mPropertyMap[ "YoungsModulusFiber" ]  = static_cast< uint >( CM_Property_Type_MT::EMOD2 ) + tCurrentIndexOffSet;
            mPropertyMap[ "PoissonRatioFiber" ]   = static_cast< uint >( CM_Property_Type_MT::NU2 ) + tCurrentIndexOffSet;
            mPropertyMap[ "VolumeFraction" ]      = static_cast< uint >( CM_Property_Type_MT::VF ) + tCurrentIndexOffSet;
            mPropertyMap[ "OrientationInPlane" ]  = static_cast< uint >( CM_Property_Type_MT::THETA_IP ) + tCurrentIndexOffSet;
            mPropertyMap[ "OrientationOutPlane" ] = static_cast< uint >( CM_Property_Type_MT::THETA_OP ) + tCurrentIndexOffSet;
            mPropertyMap[ "AspectRatio" ]         = static_cast< uint >( CM_Property_Type_MT::AR ) + tCurrentIndexOffSet;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::set_local_properties()
        {
            // set the parent properties
            CM_Struc_Linear::set_local_properties();

            // default local properties
            mPropEModFib        = this->get_property( "YoungsModulusFiber" );
            mPropPoissonFib     = this->get_property( "PoissonRatioFiber" );
            mPropEModMat        = this->get_property( "YoungsModulusMatrix" );
            mPropPoissonMat     = this->get_property( "PoissonRatioMatrix" );
            mPropVolumeFraction = this->get_property( "VolumeFraction" );
            mPropThetaIp        = this->get_property( "OrientationInPlane" );
            mPropThetaOp        = this->get_property( "OrientationOutPlane" );
            mPropAspectRatio    = this->get_property( "AspectRatio" );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::eval_const()
        {
            //  get the properties of the CM model
            const real& tNuMat   = mPropPoissonMat->val()( 0 );
            const real& tEmodMat = mPropEModMat->val()( 0 );
            const real& tEmodFib = mPropEModFib->val()( 0 );
            const real& tNuFib   = mPropPoissonFib->val()( 0 );
            const real& tAr      = mPropAspectRatio->val()( 0 );
            const real& tVf      = mPropVolumeFraction->val()( 0 );
            const real& tThetaIp = mPropThetaIp->val()( 0 );
            const real& tThetaOp = mPropThetaOp->val()( 0 );

            // evaluate the constitutive matrix
            ( this->*mConstFunc )( { tEmodMat, tEmodFib, tNuMat, tNuFib, tVf, tThetaIp, tThetaOp, tAr } );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::full_plane_stress(
                std::initializer_list< const real >&& tParams )
        {
            //  reset the constitutive tensor
            mConstPrime.fill( 0.0 );

            // get the constitutive model parameters
            const real& tEm    = tParams.begin()[ 0 ];
            const real& tEf    = tParams.begin()[ 1 ];
            const real& tNum   = tParams.begin()[ 2 ];
            const real& tNuf   = tParams.begin()[ 3 ];
            const real& tVf    = tParams.begin()[ 4 ];
            const real& theta1 = tParams.begin()[ 5 ];
            // const real& theta2 = tParams.begin()[6];
            const real& tAspectRatio = tParams.begin()[ 7 ];

            // High aspect ratio => continuous fiber
            if ( tAspectRatio > 1000.0 || tAspectRatio == 0.0 )
            {
                real tC11 =
                    ( std::pow( tEf, 2.0 ) * ( -1.0 + tVf ) * tVf * std::pow( 1.0 + tNum, 2.0 ) * ( -1.0 + 2.0 * tNum ) +
                    std::pow( tEm, 2.0 ) * ( -1.0 + tVf ) * ( 1.0 + tVf + ( -1.0 + tVf ) * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) -
                    tEm * tEf * ( 1.0 + tNum ) * ( -1.0 + tNum + tVf * ( 1.0 + tNuf - 6.0 * tNum * tNuf + tVf * ( -2.0 + tNum + tNuf + 4.0 * tNum * tNuf ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) );    // 1st row
                real tC12 =
                    ( tEm * ( -( tEf * ( 1.0 + tNum ) * ( ( -1.0 + tVf ) * tNum + 2.0 * tVf * ( -1.0 + tNum ) * tNuf ) ) +
                    tEm * ( -1.0 + tVf ) * tNum * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) ) / ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) );
                real tC22 =
                    -( ( tEm * ( -1.0 + tNum ) * ( std::pow( tEf, 2.0 ) * ( -1.0 + tVf ) * std::pow( 1.0 + tNum, 2.0 ) * ( -3.0 - 2.0 * tVf + 4.0 * ( 1.0 + tVf ) * tNum ) +
                    std::pow( tEm, 2.0 ) * ( -1.0 + tVf ) * ( 1.0 + 2.0 * tVf ) * std::pow( 1.0 + tNuf, 2.0 ) * ( -1.0 + 2.0 * tNuf ) + 2.0 * tEm * tEf * ( 1.0 + tNum ) * ( 1.0 + tNuf ) * ( 2.0 - 2.0 * tNum - 3.0 * tNuf +
                    tVf * tNuf + 4.0 * tNum * tNuf - 2.0 * std::pow( tVf, 2.0 ) * ( -1.0 + tNum + tNuf ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -3.0 + tNum + 4.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( -1.0 + tVf * ( -3.0 + 4.0 * tNum ) ) * ( 1.0 + tNuf ) ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) ) );
                real tC23 =
                    -( ( tEm * ( pow( tEf, 2.0 ) * ( -1.0 + tVf ) * pow( 1.0 + tNum, 2.0 ) * ( tVf * pow( 1.0 - 2.0 * tNum, 2.0 ) +
                    ( 3.0 - 4.0 * tNum ) * tNum ) + std::pow( tEm, 2.0 ) * ( -1.0 + tVf ) * ( -tNum + tVf * ( -1.0 + 2.0 * tNum ) ) * pow( 1.0 + tNuf, 2.0 ) * ( -1.0 + 2.0 * tNuf ) -
                    2.0 * tEm * tEf * ( 1.0 + tNum ) * ( 1.0 + tNuf ) * ( pow( tVf, 2.0 ) * ( -1.0 + 2.0 * tNum ) * ( -1.0 + tNum + tNuf ) + tVf * ( -1.0 + tNum + 5.0 * tNuf - 7.0 * tNum * tNuf ) +
                    tNum * ( 2.0 - 2.0 * tNum - 3.0 * tNuf + 4.0 * tNum * tNuf ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -3.0 + tNum + 4.0 * pow( tNum, 2.0 ) ) -
                        tEm * ( -1.0 + tVf * ( -3.0 + 4.0 * tNum ) ) * ( 1.0 + tNuf ) ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * pow( tNum, 2.0 ) ) -
                        tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * pow( tNuf, 2.0 ) ) ) ) );
                real tC33 = tC22;
                real tC66 = ( tEm * ( -( tEf * ( 1.0 + tVf ) * ( 1.0 + tNum ) ) + tEm * ( -1.0 + tVf ) * ( 1.0 + tNuf ) ) ) / ( 2.0 * ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( 1.0 + tNum ) - tEm * ( 1.0 + tVf ) * ( 1.0 + tNuf ) ) );

                // first row
                mConstPrime( 0, 0 ) = tC11 - ( tC12 * tC12 / tC33 );
                mConstPrime( 0, 1 ) = tC12 - ( tC12 * tC23 / tC33 );
                // second row
                mConstPrime( 1, 0 ) = mConstPrime( 0, 1 );
                mConstPrime( 1, 1 ) = tC22 - ( tC23 * tC23 / tC33 );
                // third row
                mConstPrime( 2, 2 ) = tC66;
            }
            else
            {
                // Eshelby Tensor
                this->set_eshelby_tensor( tAspectRatio, tNum );

                // Matrix and fiber material stiffness matrices
                this->set_isotopic_tensor( tEm, tNum, mConstMatrix );
                this->set_isotopic_tensor( tEf, tNuf, mConstFiber );

                // build an identity matrix
                Matrix< DDRMat > I;
                eye( 6, 6, I );

                // Dilute concentration factor
                // Adil = [I + S.C1^-1.(C2 - C1)]^-1
                Matrix< DDRMat > tAdil = inv( I + mEshelbyTensor * inv( mConstMatrix ) * ( mConstFiber - mConstMatrix ) );

                // Mori-Tanaka concetration tensor
                // AMT = Adil.[(1 - tVf)*I + tVf*Adil]^-1
                Matrix< DDRMat > tAMT = tAdil * inv( ( 1.0 - tVf ) * I + tVf * tAdil );

                // Get the effective stiffness matrix
                // Ceff = C1 + tVf*(C2 - C1).AMT
                Matrix< DDRMat > tCeff = mConstMatrix + tVf * ( mConstFiber - mConstMatrix ) * tAMT;

                // first row
                mConstPrime( 0, 0 ) = tCeff( 0, 0 ) - ( tCeff( 0, 1 ) * tCeff( 0, 1 ) / tCeff( 2, 2 ) );
                mConstPrime( 0, 1 ) = tCeff( 0, 1 ) - ( tCeff( 0, 1 ) * tCeff( 1, 2 ) / tCeff( 2, 2 ) );
                // second row
                mConstPrime( 1, 0 ) = mConstPrime( 0, 1 );
                mConstPrime( 1, 1 ) = tCeff( 1, 1 ) - ( tCeff( 1, 2 ) * tCeff( 1, 2 ) / tCeff( 2, 2 ) );
                // third row
                mConstPrime( 2, 2 ) = tCeff( 5, 5 );
            }

            real lx = std::cos( theta1 );
            real ly = std::sin( theta1 );

            real rx = -std::sin( theta1 );
            real ry = std::cos( theta1 );

            // first row
            mRotation( 0, 0 ) = lx * lx;
            mRotation( 0, 1 ) = ly * ly;
            mRotation( 0, 2 ) = lx * ly;
            // second row
            mRotation( 1, 0 ) = rx * rx;
            mRotation( 1, 1 ) = ry * ry;
            mRotation( 1, 2 ) = rx * ry;
            // third row
            mRotation( 2, 0 ) = 2.0 * lx * rx;
            mRotation( 2, 1 ) = 2.0 * ly * ry;
            mRotation( 2, 2 ) = lx * ry + ly * rx;

            mConst = trans( mRotation ) * mConstPrime * mRotation;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::full_3d(
                std::initializer_list< const real >&& tParams )
        {
            // reset the constitutive tensor
            mConstPrime.fill( 0.0 );

            // get the constitutive model parameters
            const real& tEm          = tParams.begin()[ 0 ];
            const real& tEf          = tParams.begin()[ 1 ];
            const real& tNum         = tParams.begin()[ 2 ];
            const real& tNuf         = tParams.begin()[ 3 ];
            const real& tVf          = tParams.begin()[ 4 ];
            const real& theta1       = tParams.begin()[ 5 ];
            const real& theta2       = tParams.begin()[ 6 ];
            const real& tAspectRatio = tParams.begin()[ 7 ];

            // High aspect ratio => continuous fiber
            if ( tAspectRatio > 1000.0 || tAspectRatio == 0.0 )
            {
                // first row
                mConstPrime( 0, 0 ) =
                    ( std::pow( tEf, 2.0 ) * ( -1.0 + tVf ) * tVf * std::pow( 1.0 + tNum, 2.0 ) * ( -1.0 + 2.0 * tNum ) +
                    std::pow( tEm, 2.0 ) * ( -1.0 + tVf ) * ( 1.0 + tVf + ( -1.0 + tVf ) * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) -
                    tEm * tEf * ( 1.0 + tNum ) * ( -1.0 + tNum + tVf * ( 1.0 + tNuf - 6.0 * tNum * tNuf + tVf * ( -2 + tNum + tNuf + 4.0 * tNum * tNuf ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) );    // 1st row
                mConstPrime( 0, 1 ) =
                    ( tEm * ( -( tEf * ( 1.0 + tNum ) * ( ( -1.0 + tVf ) * tNum + 2.0 * tVf * ( -1.0 + tNum ) * tNuf ) ) +
                    tEm * ( -1.0 + tVf ) * tNum * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) );
                mConstPrime( 0, 2 ) = mConstPrime( 0, 1 );
                // second row
                mConstPrime( 1, 0 ) = mConstPrime( 0, 1 );
                mConstPrime( 1, 1 ) =
                    -( ( tEm * ( -1.0 + tNum ) * ( std::pow( tEf, 2.0 ) * ( -1.0 + tVf ) * std::pow( 1.0 + tNum, 2.0 ) * ( -3.0 - 2.0 * tVf +
                    4.0 * ( 1.0 + tVf ) * tNum ) + std::pow( tEm, 2.0 ) * ( -1.0 + tVf ) * ( 1.0 + 2.0 * tVf ) * std::pow( 1.0 + tNuf, 2.0 ) * ( -1.0 + 2.0 * tNuf ) +
                    2.0 * tEm * tEf * ( 1.0 + tNum ) * ( 1.0 + tNuf ) * ( 2.0 - 2.0 * tNum -
                    3.0 * tNuf + tVf * tNuf + 4.0 * tNum * tNuf -
                    2.0 * std::pow( tVf, 2.0 ) * ( -1.0 + tNum + tNuf ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -3.0 + tNum + 4.0 * std::pow( tNum, 2.0 ) ) -
                    tEm * ( -1.0 + tVf * ( -3.0 + 4.0 * tNum ) ) * ( 1.0 + tNuf ) ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * std::pow( tNum, 2.0 ) ) - tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * std::pow( tNuf, 2.0 ) ) ) ) );
                mConstPrime( 1, 2 ) =
                    -( ( tEm * ( pow( tEf, 2.0 ) * ( -1.0 + tVf ) * pow( 1.0 + tNum, 2.0 ) * ( tVf * pow( 1.0 - 2.0 * tNum, 2.0 ) +
                    ( 3.0 - 4.0 * tNum ) * tNum ) + std::pow( tEm, 2.0 ) * ( -1.0 + tVf ) * ( -tNum + tVf * ( -1.0 + 2.0 * tNum ) ) * pow( 1.0 + tNuf, 2.0 ) * ( -1.0 + 2.0 * tNuf ) -
                    2.0 * tEm * tEf * ( 1.0 + tNum ) * ( 1.0 + tNuf ) * ( pow( tVf, 2.0 ) * ( -1.0 + 2.0 * tNum ) * ( -1.0 + tNum + tNuf ) +
                    tVf * ( -1.0 + tNum + 5.0 * tNuf - 7.0 * tNum * tNuf ) + tNum * ( 2.0 - 2.0 * tNum - 3.0 * tNuf + 4.0 * tNum * tNuf ) ) ) ) /
                        ( ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -3.0 + tNum + 4.0 * pow( tNum, 2.0 ) ) -
                    tEm * ( -1.0 + tVf * ( -3.0 + 4.0 * tNum ) ) * ( 1.0 + tNuf ) ) * ( tEf * ( -1.0 + tVf ) * ( -1.0 + tNum + 2.0 * pow( tNum, 2.0 ) ) -
                    tEm * ( 1.0 + tVf - 2.0 * tNum ) * ( -1.0 + tNuf + 2.0 * pow( tNuf, 2.0 ) ) ) ) );
                // third row
                mConstPrime( 2, 0 ) = mConstPrime( 0, 2 );
                mConstPrime( 2, 1 ) = mConstPrime( 1, 2 );
                mConstPrime( 2, 2 ) = mConstPrime( 1, 1 );
                // shear terms - order  xy, xz, yz (older version had yz, xz, xy)
                mConstPrime( 5, 5 ) =
                    ( tEm * ( tEf * ( 3.0 + tVf - 4.0 * tNum ) * ( 1.0 + tNum ) -
                    tEm * ( -1.0 + tVf ) * ( 1.0 + tNuf ) ) ) / ( 2.0 * ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( -3.0 + tNum + 4.0 * pow( tNum, 2.0 ) ) -
                    tEm * ( -1.0 + tVf * ( -3.0 + 4.0 * tNum ) ) * ( 1.0 + tNuf ) ) );
                mConstPrime( 4, 4 ) =
                    ( tEm * ( -( tEf * ( 1.0 + tVf ) * ( 1.0 + tNum ) ) +
                    tEm * ( -1.0 + tVf ) * ( 1.0 + tNuf ) ) ) / ( 2.0 * ( 1.0 + tNum ) * ( tEf * ( -1.0 + tVf ) * ( 1.0 + tNum ) -
                    tEm * ( 1.0 + tVf ) * ( 1.0 + tNuf ) ) );
                mConstPrime( 3, 3 ) = mConstPrime( 4, 4 );
            }
            else    // Low aspect ratio, particles
            {
                // Eshelby Tensor
                this->set_eshelby_tensor( tAspectRatio, tNum );

                // Matrix and fiber material stiffness matrices
                this->set_isotopic_tensor( tEm, tNum, mConstMatrix );
                this->set_isotopic_tensor( tEf, tNuf, mConstFiber );

                // build an identity matrix
                Matrix< DDRMat > I;
                eye( 6, 6, I );

                // Dilute concentration factor
                // Adil = [I + S.C1^-1.(C2 - C1)]^-1
                Matrix< DDRMat > tAdil = inv( I + mEshelbyTensor * inv( mConstMatrix ) * ( mConstFiber - mConstMatrix ) );

                // Mori-Tanaka concentration tensor
                // AMT = Adil.[(1 - tVf)*I + tVf*Adil]^-1
                Matrix< DDRMat > tAMT = tAdil * inv( ( 1.0 - tVf ) * I + tVf * tAdil );

                // Get the effective stiffness matrix
                // CeffMT = C1 + tVf*(C2 - C1).AMT
                mConstPrime = mConstMatrix + tVf * ( mConstFiber - mConstMatrix ) * tAMT;

                // Exchange C44 with C66
                std::swap( mConstPrime( 3, 3 ), mConstPrime( 5, 5 ) );
            }

            real lx = std::cos( theta2 ) * std::cos( theta1 );
            real ly = std::sin( theta1 );
            real lz = std::sin( theta2 ) * std::cos( theta1 );

            real rx = -std::cos( theta2 ) * std::sin( theta1 );
            real ry = std::cos( theta1 );
            real rz = -std::sin( theta2 ) * std::sin( theta1 );

            real tx = -std::sin( theta2 );
            real ty = 0.0;
            real tz = std::cos( theta2 );

            // first row
            mRotation( 0, 0 ) = lx * lx;
            mRotation( 0, 1 ) = ly * ly;
            mRotation( 0, 2 ) = lz * lz;
            mRotation( 0, 3 ) = ly * lz;
            mRotation( 0, 4 ) = lx * lz;
            mRotation( 0, 5 ) = lx * ly;
            // second row
            mRotation( 1, 0 ) = rx * rx;
            mRotation( 1, 1 ) = ry * ry;
            mRotation( 1, 2 ) = rz * rz;
            mRotation( 1, 3 ) = ry * rz;
            mRotation( 1, 4 ) = rx * rz;
            mRotation( 1, 5 ) = rx * ry;
            // third row
            mRotation( 2, 0 ) = tx * tx;
            mRotation( 2, 1 ) = ty * ty;
            mRotation( 2, 2 ) = tz * tz;
            mRotation( 2, 3 ) = ty * tz;
            mRotation( 2, 4 ) = tx * tz;
            mRotation( 2, 5 ) = tx * ty;
            // fourth row
            mRotation( 3, 0 ) = 2.0 * lx * rx;
            mRotation( 3, 1 ) = 2.0 * ly * ry;
            mRotation( 3, 2 ) = 2.0 * lz * rz;
            mRotation( 3, 3 ) = ly * rz + lz * ry;
            mRotation( 3, 4 ) = lz * rx + lx * rz;
            mRotation( 3, 5 ) = lx * ry + ly * rx;
            // fifth row
            mRotation( 4, 0 ) = 2.0 * lx * tx;
            mRotation( 4, 1 ) = 2.0 * ly * ty;
            mRotation( 4, 2 ) = 2.0 * lz * tz;
            mRotation( 4, 3 ) = ty * lz + tz * ly;
            mRotation( 4, 4 ) = tz * lx + tx * lz;
            mRotation( 4, 5 ) = tx * ly + ty * lx;
            // sixth row
            mRotation( 5, 0 ) = 2.0 * rx * tx;
            mRotation( 5, 1 ) = 2.0 * ry * ty;
            mRotation( 5, 2 ) = 2.0 * rz * tz;
            mRotation( 5, 3 ) = ry * tz + rz * ty;
            mRotation( 5, 4 ) = rz * tx + rx * tz;
            mRotation( 5, 5 ) = rx * ty + ry * tx;

            mConst = trans( mRotation ) * mConstPrime * mRotation;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::set_eshelby_tensor( real const & tAspectRatio, real const & v )
        {
            moris::real tAspectRatio2 = tAspectRatio * tAspectRatio;

            real g = 1.0;
            if ( tAspectRatio == 1.0 )
            {
                // Eshelby Tensor
                mEshelbyTensor( 0, 0 ) = ( 7.0 - 5.0 * v ) / ( 15.0 - 15.0 * v );
                mEshelbyTensor( 0, 1 ) = ( 1.0 - 5.0 * v ) / ( 15.0 * ( -1.0 + v ) );
                mEshelbyTensor( 0, 2 ) = mEshelbyTensor( 0, 1 );

                mEshelbyTensor( 1, 0 ) = mEshelbyTensor( 0, 1 );
                mEshelbyTensor( 1, 1 ) = mEshelbyTensor( 0, 0 );
                mEshelbyTensor( 1, 2 ) = mEshelbyTensor( 0, 1 );

                mEshelbyTensor( 2, 0 ) = mEshelbyTensor( 0, 1 );
                mEshelbyTensor( 2, 1 ) = mEshelbyTensor( 2, 0 );
                mEshelbyTensor( 2, 2 ) = mEshelbyTensor( 0, 0 );

                mEshelbyTensor( 3, 3 ) = ( 2.0 * ( -4.0 + 5.0 * v ) ) / ( 15.0 * ( -1.0 + v ) );
                mEshelbyTensor( 4, 4 ) = mEshelbyTensor( 3, 3 );
                mEshelbyTensor( 5, 5 ) = mEshelbyTensor( 3, 3 );

                return;
            }

            if ( tAspectRatio < 1.0 )
            {
                g = ( tAspectRatio * ( std::acos( tAspectRatio ) - tAspectRatio * std::pow( 1.0 - tAspectRatio2, 0.5 ) ) ) /
                    std::pow( 1.0 - tAspectRatio2, 1.5 );
            }

            if ( tAspectRatio > 1.0 )
            {
                g = ( tAspectRatio * ( tAspectRatio * std::pow( tAspectRatio2 - 1, 0.5 ) - std::acosh( tAspectRatio ) ) ) /
                    std::pow( tAspectRatio2 - 1.0, 1.5 );
            }

            // Eshelby Tensor
            mEshelbyTensor( 0, 0 ) =
                ( 1.0 + ( -1.0 + 3.0 * tAspectRatio2 ) / ( -1.0 + tAspectRatio2 ) - 2.0 * v ) / ( 2.0 * ( 1.0 - v ) ) -
                ( ( 1.0 + ( 3.0 * tAspectRatio2 ) / ( -1.0 + tAspectRatio2 ) - 2.0 * v ) * g ) / ( 2.0 * ( 1.0 - v ) );
            mEshelbyTensor( 0, 1 ) =
                -( 1.0 + 1.0 / ( -1.0 + tAspectRatio2 ) - 2.0 * v ) / ( 2.0 * ( 1.0 - v ) ) +
                ( ( 1.0 + 3.0 / ( 2.0 * ( -1.0 + tAspectRatio2 ) ) - 2.0 * v ) * g ) / ( 2.0 * ( 1.0 - v ) );
            mEshelbyTensor( 0, 2 ) = mEshelbyTensor( 0, 1 );

            mEshelbyTensor( 1, 0 ) =
                -tAspectRatio2 / ( 2.0 * ( -1.0 + tAspectRatio2 ) * ( 1.0 - v ) ) +
                ( ( -1.0 + ( 3.0 * tAspectRatio2 ) / ( -1.0 + tAspectRatio2 ) + 2.0 * v ) * g ) / ( 4.0 * ( 1.0 - v ) );
            mEshelbyTensor( 1, 1 ) =
                ( 3.0 * tAspectRatio2 ) / ( 8.0 * ( -1.0 + tAspectRatio2 ) * ( 1.0 - v ) ) +
                ( ( 1.0 - 9.0 / ( 4.0 * ( -1.0 + tAspectRatio2 ) ) - 2.0 * v ) * g ) / ( 4.0 * ( 1.0 - v ) );
            mEshelbyTensor( 1, 2 ) =
                tAspectRatio2 / ( 8.0 * ( -1.0 + tAspectRatio2 ) * ( 1.0 - v ) ) -
                ( ( 1.0 + 3.0 / ( 4.0 * ( -1.0 + tAspectRatio2 ) ) - 2.0 * v ) * g ) / ( 4.0 * ( 1.0 - v ) );

            mEshelbyTensor( 2, 0 ) = mEshelbyTensor( 1, 0 );
            mEshelbyTensor( 2, 1 ) = mEshelbyTensor( 1, 2 );
            mEshelbyTensor( 2, 2 ) = mEshelbyTensor( 1, 1 );

            mEshelbyTensor( 3, 3 ) =
                tAspectRatio2 / ( 4.0 * ( -1.0 + tAspectRatio2 ) * ( 1.0 - v ) ) +
                ( ( 1.0 - 3.0 / ( 4.0 * ( -1.0 + tAspectRatio2 ) ) - 2.0 * v ) * g ) / ( 2.0 * ( 1.0 - v ) );
            mEshelbyTensor( 4, 4 ) =
                ( 1.0 - ( 1.0 + tAspectRatio2 ) / ( -1.0 + tAspectRatio2 ) - 2.0 * v ) / ( 2.0 * ( 1.0 - v ) ) -
                ( ( 1.0 - ( 3.0 * ( 1.0 + tAspectRatio2 ) ) / ( -1.0 + tAspectRatio2 ) - 2.0 * v ) * g ) / ( 4.0 * ( 1.0 - v ) );
            mEshelbyTensor( 5, 5 ) = mEshelbyTensor( 4, 4 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::set_isotopic_tensor( const real& aEmod,
                const real&                                          aNu,
                Matrix< DDRMat >&                                    aIsotropicConst )
        {
            const real tPre = aEmod / ( 1.0 + aNu ) / ( 1.0 - 2.0 * aNu );

            aIsotropicConst( 0, 0 ) = tPre * ( 1.0 - aNu );
            aIsotropicConst( 0, 1 ) = tPre * aNu;
            aIsotropicConst( 0, 2 ) = tPre * aNu;
            aIsotropicConst( 1, 0 ) = tPre * aNu;
            aIsotropicConst( 1, 1 ) = tPre * ( 1.0 - aNu );
            aIsotropicConst( 1, 2 ) = tPre * aNu;
            aIsotropicConst( 2, 0 ) = tPre * aNu;
            aIsotropicConst( 2, 1 ) = tPre * aNu;
            aIsotropicConst( 2, 2 ) = tPre * ( 1.0 - aNu );
            aIsotropicConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
            aIsotropicConst( 4, 4 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
            aIsotropicConst( 5, 5 ) = tPre * ( 1.0 - 2.0 * aNu ) / 2.0;
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::set_function_pointers()
        {
            CM_Struc_Linear::set_function_pointers();

            switch ( mSpaceDim )
            {
                case 2:
                {
                    switch ( mPlaneType )
                    {
                        case Model_Type::PLANE_STRESS:
                        {
                            mRotation.set_size( 3, 3, 0.0 );
                            mConstPrime.set_size( 3, 3, 0.0 );
                            mRotationDerInPlane.set_size( 3, 3, 0.0 );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false,
                                    "Mori Tanaka in 2d requires "
                                    "plane stress" );
                        }
                    }
                    break;
                }

                case 3:
                {
                    mRotation.set_size( 6, 6, 0.0 );
                    mConstPrime.set_size( 6, 6, 0.0 );
                    mRotationDerInPlane.set_size( 6, 6, 0.0 );
                    mRotationDerOutPlane.set_size( 6, 6, 0.0 );
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "Mori Tanaka implemented for 2d and 3d only" );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
        {
            // call the parent contribution
            CM_Struc_Linear::eval_dFluxdDOF( aDofTypes );

            // get the dof type as a uint
            const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // if elastic modulus depends on dof type
            if ( mPropThetaIp->check_dof_dependency( aDofTypes ) )
            {
                // update constitutive matrix and rotation tensor
                this->eval_const();

                // evaluate the derivative of the rotation tensor
                this->eval_inplane_rotation_derivative();

                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ) +=
                        ( trans( mRotationDerInPlane ) * mConstPrime * mRotation              //
                                + trans( mRotation ) * mConstPrime * mRotationDerInPlane )    //
                        * this->strain() * mPropThetaIp->dPropdDOF( aDofTypes );
            }

            // if elastic modulus depends on dof type
            if ( mPropThetaOp->check_dof_dependency( aDofTypes ) && mSpaceDim == 3 )
            {
                // update constitutive matrix and rotation tensor
                this->eval_const();

                // evaluate the derivative of the rotation tensor
                this->eval_outplane_rotation_derivative();

                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ) +=
                        ( trans( mRotationDerOutPlane ) * mConstPrime * mRotation              //
                                + trans( mRotation ) * mConstPrime * mRotationDerOutPlane )    //
                        * this->strain() * mPropThetaOp->dPropdDOF( aDofTypes );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::eval_dTestTractiondDOF(
                const Vector< MSI::Dof_Type >& aDofTypes,
                const Matrix< DDRMat >&      aNormal,
                const Matrix< DDRMat >&      aJump,
                const Vector< MSI::Dof_Type >& aTestDofTypes )
        {
            CM_Struc_Linear::eval_dTestTractiondDOF( aDofTypes, aNormal, aJump, aTestDofTypes );

            // get test dof type index
            const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof type index
            const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if test traction wrt displacement
            if ( aTestDofTypes( 0 ) == mDofDispl )
            {
                // if elastic modulus depends on dof type
                if ( mPropThetaIp->check_dof_dependency( aDofTypes ) )
                {
                    // update consitutive matrix and rotation tensor
                    this->eval_const();

                    // evalute the derivative of the rotation tensor
                    this->eval_inplane_rotation_derivative();

                    // flatten the normal
                    Matrix< DDRMat > tFlatNormal;
                    this->flatten_normal( aNormal, tFlatNormal );

                    // compute test traction wrt displacement
                    mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                            trans( this->testStrain() )
                            * ( trans( mRotationDerInPlane ) * mConstPrime * mRotation + trans( mRotation ) * mConstPrime * mRotationDerInPlane )
                            * trans( tFlatNormal ) * aJump * mPropThetaIp->dPropdDOF( aDofTypes );
                }

                // if elastic modulus depends on dof type
                if ( mPropThetaOp->check_dof_dependency( aDofTypes ) && mSpaceDim == 3 )
                {
                    // update consitutive matrix and rotation tensor
                    this->eval_const();

                    // evalute the derivative of the rotation tensor
                    this->eval_outplane_rotation_derivative();

                    // flatten the normal
                    Matrix< DDRMat > tFlatNormal;
                    this->flatten_normal( aNormal, tFlatNormal );

                    // compute test traction wrt displacement
                    mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                            trans( this->testStrain() )
                            * ( trans( mRotationDerOutPlane ) * mConstPrime * mRotation + trans( mRotation ) * mConstPrime * mRotationDerOutPlane )
                            * trans( tFlatNormal ) * aJump * mPropThetaOp->dPropdDOF( aDofTypes );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::eval_inplane_rotation_derivative()
        {
            switch ( mSpaceDim )
            {
                case 2:
                {
                    // get rotation angle
                    const real& theta1 = mPropThetaIp->val()( 0 );

                    // compute derivative of ration matrix wrt in-plane angle
                    real lx2 = std::cos( 2.0 * theta1 );
                    real ly2 = std::sin( 2.0 * theta1 );

                    // first row
                    mRotationDerInPlane( 0, 0 ) = -ly2;
                    mRotationDerInPlane( 0, 1 ) = ly2;
                    mRotationDerInPlane( 0, 2 ) = lx2;
                    // second row
                    mRotationDerInPlane( 1, 0 ) = ly2;
                    mRotationDerInPlane( 1, 1 ) = -ly2;
                    mRotationDerInPlane( 1, 2 ) = -lx2;
                    // third row
                    mRotationDerInPlane( 2, 0 ) = -2.0 * lx2;
                    mRotationDerInPlane( 2, 1 ) = 2.0 * lx2;
                    mRotationDerInPlane( 2, 2 ) = -2.0 * ly2;

                    break;
                }

                case 3:
                {
                    // get rotation angle
                    const real& thetai = mPropThetaIp->val()( 0 );
                    const real& thetao = mPropThetaOp->val()( 0 );

                    // name cos and sin for convience
                    real ci  = std::cos( thetai );
                    real si  = std::sin( thetai );
                    real co  = std::cos( thetao );
                    real so  = std::sin( thetao );
                    real c2i = std::cos( 2.0 * thetai );
                    real s2i = std::sin( 2.0 * thetai );
                    real c2o = std::cos( 2.0 * thetao );
                    real s2o = std::sin( 2.0 * thetao );

                    // first row
                    mRotationDerInPlane( 0, 0 ) = -2.0 * ci * co * co * si;
                    mRotationDerInPlane( 0, 1 ) = s2i;
                    mRotationDerInPlane( 0, 2 ) = -2.0 * ci * si * so * so;
                    mRotationDerInPlane( 0, 3 ) = c2i * so;
                    mRotationDerInPlane( 0, 4 ) = -2.0 * ci * co * si * so;
                    mRotationDerInPlane( 0, 5 ) = c2i * co;

                    // second row
                    mRotationDerInPlane( 1, 0 ) = 2.0 * ci * co * co * si;
                    mRotationDerInPlane( 1, 1 ) = -s2i;
                    mRotationDerInPlane( 1, 2 ) = 2.0 * ci * si * so * so;
                    mRotationDerInPlane( 1, 3 ) = -c2i * so;
                    mRotationDerInPlane( 1, 4 ) = 2.0 * ci * co * si * so;
                    mRotationDerInPlane( 1, 5 ) = -c2i * co;

                    // third row
                    mRotationDerInPlane( 2, 0 ) = 0.0;
                    mRotationDerInPlane( 2, 1 ) = 0.0;
                    mRotationDerInPlane( 2, 2 ) = 0.0;
                    mRotationDerInPlane( 2, 3 ) = 0.0;
                    mRotationDerInPlane( 2, 4 ) = 0.0;
                    mRotationDerInPlane( 2, 5 ) = 0.0;

                    // fourth row
                    mRotationDerInPlane( 3, 0 ) = -2.0 * c2i * co * co;
                    mRotationDerInPlane( 3, 1 ) = 2.0 - 4.0 * si * si;
                    mRotationDerInPlane( 3, 2 ) = 2.0 * so * so * ( 2.0 * si * si - 1 );
                    mRotationDerInPlane( 3, 3 ) = -4.0 * ci * si * so;
                    mRotationDerInPlane( 3, 4 ) = -c2i * s2o;
                    mRotationDerInPlane( 3, 5 ) = -4.0 * ci * co * si;

                    // fifth row
                    mRotationDerInPlane( 4, 0 ) = 2.0 * co * si * so;
                    mRotationDerInPlane( 4, 1 ) = 0.0;
                    mRotationDerInPlane( 4, 2 ) = -2.0 * co * si * so;
                    mRotationDerInPlane( 4, 3 ) = ci * co;
                    mRotationDerInPlane( 4, 4 ) = -c2o * si;
                    mRotationDerInPlane( 4, 5 ) = -ci * so;

                    // sixth row
                    mRotationDerInPlane( 5, 0 ) = 2.0 * ci * co * so;
                    mRotationDerInPlane( 5, 1 ) = 0.0;
                    mRotationDerInPlane( 5, 2 ) = -2.0 * ci * co * so;
                    mRotationDerInPlane( 5, 3 ) = -co * si;
                    mRotationDerInPlane( 5, 4 ) = -c2o * ci;
                    mRotationDerInPlane( 5, 5 ) = si * so;

                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "Mori Tanaka implemented for 2d and 3d only" );
                    break;
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::eval_outplane_rotation_derivative()
        {
            // get rotation angle
            const real& thetai = mPropThetaIp->val()( 0 );
            const real& thetao = mPropThetaOp->val()( 0 );

            // name cos and sin for convenience
            real ci  = std::cos( thetai );
            real si  = std::sin( thetai );
            real co  = std::cos( thetao );
            real so  = std::sin( thetao );
            real c2i = std::cos( 2.0 * thetai );
            real s2i = std::sin( 2.0 * thetai );
            real c2o = std::cos( 2.0 * thetao );
            real s2o = std::sin( 2.0 * thetao );

            // first row
            mRotationDerOutPlane( 0, 0 ) = -2.0 * ci * ci * co * so;
            mRotationDerOutPlane( 0, 1 ) = 0.0;
            mRotationDerOutPlane( 0, 2 ) = 2.0 * ci * ci * co * so;
            mRotationDerOutPlane( 0, 3 ) = ci * co * si;
            mRotationDerOutPlane( 0, 4 ) = c2o * ci * ci;
            mRotationDerOutPlane( 0, 5 ) = -ci * si * so;

            // second row
            mRotationDerOutPlane( 1, 0 ) = -2.0 * co * si * si * so;
            mRotationDerOutPlane( 1, 1 ) = 0.0;
            mRotationDerOutPlane( 1, 2 ) = 2.0 * co * si * si * so;
            mRotationDerOutPlane( 1, 3 ) = -ci * co * si;
            mRotationDerOutPlane( 1, 4 ) = -si * si * ( 2.0 * so * so - 1.0 );
            mRotationDerOutPlane( 1, 5 ) = ci * si * so;

            // third row
            mRotationDerOutPlane( 2, 0 ) = s2o;
            mRotationDerOutPlane( 2, 1 ) = 0.0;
            mRotationDerOutPlane( 2, 2 ) = -s2o;
            mRotationDerOutPlane( 2, 3 ) = 0.0;
            mRotationDerOutPlane( 2, 4 ) = 2.0 * so * so - 1;
            mRotationDerOutPlane( 2, 5 ) = 0.0;

            // fourth row
            mRotationDerOutPlane( 3, 0 ) = 4.0 * ci * co * si * so;
            mRotationDerOutPlane( 3, 1 ) = 0.0;
            mRotationDerOutPlane( 3, 2 ) = -4.0 * ci * co * si * so;
            mRotationDerOutPlane( 3, 3 ) = c2i * co;
            mRotationDerOutPlane( 3, 4 ) = -c2o * s2i;
            mRotationDerOutPlane( 3, 5 ) = -c2i * so;

            // fifth row
            mRotationDerOutPlane( 4, 0 ) = -2.0 * c2o * ci;
            mRotationDerOutPlane( 4, 1 ) = 0.0;
            mRotationDerOutPlane( 4, 2 ) = 2.0 * c2o * ci;
            mRotationDerOutPlane( 4, 3 ) = -si * so;
            mRotationDerOutPlane( 4, 4 ) = -4.0 * ci * co * so;
            mRotationDerOutPlane( 4, 5 ) = -co * si;

            // sixth row
            mRotationDerOutPlane( 5, 0 ) = 2.0 * c2o * si;
            mRotationDerOutPlane( 5, 1 ) = 0.0;
            mRotationDerOutPlane( 5, 2 ) = -2.0 * c2o * si;
            mRotationDerOutPlane( 5, 3 ) = -ci * so;
            mRotationDerOutPlane( 5, 4 ) = 4.0 * co * si * so;
            mRotationDerOutPlane( 5, 5 ) = -ci * co;
        }
        //--------------------------------------------------------------------------------------------------------------

        void
        CM_Struc_Linear_MoriTanaka::eval_dConstdDOF(
                const moris::Vector< MSI::Dof_Type >& aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator* tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdConstdDof( tDofIndex ).set_size( mConst.n_rows(), tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // if conductivity depends on the dof type
            if ( mPropThetaIp->check_dof_dependency( aDofTypes ) )
            {
                // update consitutive matrix and rotation tensor
                this->eval_const();

                // evalute the derivative of the rotation tensor
                this->eval_inplane_rotation_derivative();

                Matrix< DDRMat > tdConstdProp = trans( mRotationDerInPlane ) * mConstPrime * mRotation    //
                                              + trans( mRotation ) * mConstPrime * mRotationDerInPlane;

                // compute derivative with indirect dependency through properties
                mdConstdDof( tDofIndex ) = tdConstdProp.get_column( 0 )    //
                                         * mPropThetaIp->dPropdDOF( aDofTypes );
            }
        }
    }    // namespace fem
}    // namespace moris

