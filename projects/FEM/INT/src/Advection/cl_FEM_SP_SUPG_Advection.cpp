
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
            std::shared_ptr< Property > & tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // get the norm of velocity
            real tNorm = norm( tVelocityFI->val() );

            // get the abs term
            real tAbs = 0.0;
            uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tAdd =
                        trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                tAbs += std::abs( tAdd( 0, 0 ) );
            }

            // compute hugn
            real tHugn = 1.0;
            if( tAbs > 0.0 && tNorm > 0.0 )
            {
                tHugn = 2.0 * tNorm / tAbs;
            }

            // compute tau1
            real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            real tTau2 = 4.0 * tPropConductivity->val()( 0 ) / std::pow( tHugn, 2.0 );

            // compute time increment deltat
            real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // compute tau3
            real tTau3 = 2 / tDeltaT;

            // compute stabilization parameter value
            mPPVal = {{ std::pow(
                    std::pow( tTau1, 2.0 ) +
                    std::pow( tTau2, 2.0 ) +
                    std::pow( tTau3, 2.0 ) , -0.5 ) }};

        }

        //------------------------------------------------------------------------------

        void SP_SUPG_Advection::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof, dTau1dDof, dTau3dDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
            Matrix< DDRMat > tdTau1dDof( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
            Matrix< DDRMat > tdTau2dDof( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the conductivity property
            std::shared_ptr< Property > & tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );

            // compute the norm of velocity
            real tNorm = norm( tVelocityFI->val() );

            // compute the abs term
            real tAbs = 0.0;
            uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tAdd = trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                tAbs += std::abs( tAdd( 0, 0 ) );
            }

            // compute hugn
            real tHugn = 1.0;
            if( tAbs > 0.0 && tNorm > 0.0 )
            {
                tHugn = 2.0 * tNorm / tAbs;
            }

            // compute tau1
            real tTau1 = 2.0 * tNorm / tHugn;

            // compute tau2
            real tTau2 = 4.0 * tPropConductivity->val()( 0 ) / std::pow( tHugn, 2.0 );

            // if dof type is velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                // if tAbs and tNorm > 0.0
                if( tAbs > 0.0 && tNorm > 0.0 )
                {
                    // compute derivative of hugn wrt velocity dof
                    Matrix< DDRMat > tdHugndu( 1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );
                    Matrix< DDRMat > tdNormdu( 1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );
                    Matrix< DDRMat > tdAbsdu(  1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );

                    // compute derivative of the velocity norm
                    tdNormdu +=
                            trans( tVelocityFI->val() ) * tVelocityFI->N() / tNorm;

                    // compute derivative of the abs term
                    uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
                    for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                    {
                        Matrix< DDRMat > tAdd =
                                trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
                        tdAbsdu +=
                                tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) *
                                tVelocityFI->N() / std::abs( tAdd( 0, 0 ) );
                    }

                    // compute derivative of hugn
                    tdHugndu +=
                            2.0 * ( tdNormdu * tAbs - tdAbsdu * tNorm ) / std::pow( tAbs, 2.0 );

                    // compute dtau1du
                    tdTau1dDof +=
                            2.0 * ( tHugn * tdNormdu - tdHugndu * tNorm ) / std::pow( tHugn, 2 );

                    // compute dtau2du
                    tdTau2dDof -=
                            8.0 * tPropConductivity->val()( 0 ) * tdHugndu / std::pow( tHugn, 3 );
                }
            }

            // if conductivity property depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute dtau3du
                tdTau2dDof +=
                        4.0 * tPropConductivity->dPropdDOF( aDofTypes ) / std::pow( tHugn, 2.0 );
            }

            // evaluate tau
            Matrix< DDRMat > tTau = this->val();

            // scale dSPdu
            mdPPdMasterDof( tDofIndex ) -= 0.5 * std::pow( tTau( 0 ), 3 ) *
                    ( 2.0 * tTau1 * tdTau1dDof + 2.0 * tTau2 * tdTau2dDof );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


