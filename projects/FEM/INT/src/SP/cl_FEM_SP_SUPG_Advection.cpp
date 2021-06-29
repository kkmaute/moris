
//FEM/INT/src
#include "cl_FEM_SP_SUPG_Advection.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        SP_SUPG_Advection::SP_SUPG_Advection()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Conductivity" ] = static_cast< uint >( Property_Type::CONDUCTIVITY );
            mPropertyMap[ "Source" ]       = static_cast< uint >( Property_Type::SOURCE );

        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // set dof type list
                    mMasterDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( tDofString == "Velocity" )
                        {
                            mMasterDofVelocity = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_SUPG_Advection::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_SUPG_Advection::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the conductivity property
            const std::shared_ptr< Property > & tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the mass source
            const std::shared_ptr< Property > & tSourceProp =
                    mMasterProp( static_cast< uint >( Property_Type::SOURCE ) );

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon);

            // get the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real tAbs = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tAdd =
                        trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                tAbs += std::abs( tAdd( 0, 0 ) );
            }

            // threshold tAbs
            tAbs = std::max(tAbs, mEpsilon);

            // compute and threshold hugn
            const real tHugn = std::max( 2.0 * tNorm / tAbs, mEpsilon);

            // compute tau1
            const real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            const real tTau2 = 4.0 * tPropConductivity->val()( 0 ) / std::pow( tHugn, 2.0 );

            // compute time increment deltat
            const real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute tau3
            const real tTau3 = 2.0 / tDeltaT;

            // compute sum of square terms
            real tSum = std::pow( tTau1, 2.0 ) + std::pow( tTau2, 2.0 ) + std::pow( tTau3, 2.0 );

            // add contribution from source term
            if ( tSourceProp != nullptr )
            {
                tSum += std::pow( tSourceProp->val()(0), 2.0 );
            }

            // threshold sum of square terms
            tSum = std::max(tSum, mEpsilon);

            // compute stabilization parameter value
            mPPVal = {{ std::pow( tSum, -0.5 ) }};
        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            const uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof, dTau1dDof, dTau3dDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            Matrix< DDRMat > tdTau1dDof( 1, tFIDer->get_number_of_space_time_coefficients() );
            Matrix< DDRMat > tdTau2dDof( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the conductivity property
            const std::shared_ptr< Property > & tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the mass source
            const std::shared_ptr< Property > & tSourceProp =
                    mMasterProp( static_cast< uint >( Property_Type::SOURCE ) );

            // compute and threshold the velocity norm (thresholding for consistency with derivatives)
            const real tNorm = std::max( norm( tVelocityFI->val() ), mEpsilon);

            // compute the abs term
            const uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            real tAbs = 0.0;

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tAdd = trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                tAbs += std::abs( tAdd( 0, 0 ) );
            }

            // threshold tAbs
            tAbs = std::max(tAbs, mEpsilon);

            // compute and threshold hugn
            const real tHugn = std::max( 2.0 * tNorm / tAbs, mEpsilon);

            // compute tau1
            const real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            const real tTau2 = 4.0 * tPropConductivity->val()( 0 ) / std::pow( tHugn, 2.0 );

            // compute time increment deltat
            const real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute tau3
            const real tTau3 = 2.0 / tDeltaT;

            // if dof type is velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                // compute derivative of hugn wrt velocity dof
                Matrix< DDRMat > tdHugndu( 1, tVelocityFI->get_number_of_space_time_coefficients() );
                Matrix< DDRMat > tdNormdu( 1, tVelocityFI->get_number_of_space_time_coefficients() );
                Matrix< DDRMat > tdAbsdu(  1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );

                // compute derivative of the velocity norm (compute only derivative if not thresholded)
                if ( tNorm > mEpsilon )
                {
                    tdNormdu = trans( tVelocityFI->val() ) * tVelocityFI->N() / tNorm;
                }
                else
                {
                    tdNormdu.fill( 0.0 );
                }

                // compute derivative of the abs term (compute only derivative if not thresholded)
                if ( tAbs > mEpsilon )
                {
                    uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();

                    for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                    {
                        const Matrix< DDRMat > tAdd =
                                trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );

                        // handle case that tAdd( 0, 0 ) is smaller than threshold
                        if ( std::abs( tAdd( 0, 0 ) ) > mEpsilon )
                        {
                            tdAbsdu +=
                                    tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) *
                                    tVelocityFI->N() / std::abs( tAdd( 0, 0 ) );
                        }
                    }
                }

                // compute derivative of hugn (compute only derivative if not thresholded)
                if ( tHugn > mEpsilon )
                {
                    tdHugndu = 2.0 * ( tdNormdu * tAbs - tdAbsdu * tNorm ) / std::pow( tAbs, 2.0 );
                }
                else
                {
                    tdHugndu.fill( 0.0 );
                }

                // compute dtau1du
                tdTau1dDof = 2.0 * ( tHugn * tdNormdu - tdHugndu * tNorm ) / std::pow( tHugn, 2.0 );

                // compute dtau2du
                tdTau2dDof = - 8.0 * tPropConductivity->val()( 0 ) * tdHugndu / std::pow( tHugn, 3.0 );
            }
            else
            {
                tdTau1dDof.fill( 0.0 );
                tdTau2dDof.fill( 0.0 );
            }

            // if conductivity property depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute dtau3du
                tdTau2dDof += 4.0 * tPropConductivity->dPropdDOF( aDofTypes ) / std::pow( tHugn, 2.0 );
            }

            // compute sum of square terms
            real tSum = std::pow( tTau1, 2.0 ) + std::pow( tTau2, 2.0 ) + std::pow( tTau3, 2.0 );

            // add contribution from source term
            if ( tSourceProp != nullptr )
            {
                tSum += std::pow( tSourceProp->val()(0), 2.0 );
            }

            // compute dSPdu
            if ( tSum > mEpsilon )
            {
                const real tPrefactor = - std::pow( tSum, -1.5 );

                mdPPdMasterDof( tDofIndex ) = tPrefactor * ( tTau1 * tdTau1dDof + tTau2 * tdTau2dDof);

                if ( tSourceProp != nullptr )
                {
                    if( tSourceProp->check_dof_dependency( aDofTypes ) )
                    {
                        // compute dtau3du
                        mdPPdMasterDof( tDofIndex ) +=  tPrefactor * tSourceProp->val()(0) * tSourceProp->dPropdDOF( aDofTypes );
                    }
                }
            }
            else
            {
                mdPPdMasterDof( tDofIndex ).fill( 0.0 );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mdPPdMasterDof( tDofIndex ) ),
                    "SP_SUPG_Advection::eval_dSPdMasterDOF - mdPPdMasterDof contains NAN or INF, exiting for tDofIndex = %d !\n",
                    tDofIndex);
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
