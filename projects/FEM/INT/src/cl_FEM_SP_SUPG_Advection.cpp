
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
            mPropertyMap[ "Conductivity" ]     = Property_Type::CONDUCTIVITY;
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
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_SUPG_Advection::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
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
                    break;
            }
        }


        //------------------------------------------------------------------------------
        void SP_SUPG_Advection::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_SUPG_Advection::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_SUPG_Advection::eval_SP()
        {
            // set size for SP values
            mPPVal.set_size( 1, 1, 0.0 );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // compute norm( v )
            real tNormA = std::sqrt( dot( tVelocityFI->val(), tVelocityFI->val() ) );

            // compute the diffusion coefficient
            std::shared_ptr< Property > tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );
            real tK = tPropConductivity->val()( 0 );

            // evaluate tau
            real tTau =
                    std::pow( 2.0 * tNormA / mElementSize, 2 ) +
                    std::pow( 4.0 * tK / std::pow( mElementSize, 2.0 ), 2 );

            // set tau
            mPPVal = {{ std::pow( tTau, -0.5 ) }};

//            // get the velocity FI
//            Field_Interpolator * tVelocityFI =
//                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );
//
//            // get the conductivity property
//            std::shared_ptr< Property > tPropConductivity =
//                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );
//
//            // get the norm of velocity
//            real tNorm = norm( tVelocityFI->val() );
//
//            // get the abs term
//            real tAbs = 0.0;
//            uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
//            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
//            {
//                Matrix< DDRMat > tAdd =
//                        trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
//                tAbs += std::abs( tAdd( 0, 0 ) );
//            }
//
//            // get hugn
//            real tHugn = 2.0 * tNorm / tAbs;
//
//            // get tau1
//            real tTau1 = 0.5 * tHugn / tNorm;
//
//            // compute time increment deltat
//            Matrix< DDRMat > tTimeCoeff =
//                    mMasterFIManager->get_IP_geometry_interpolator()->get_time_coeff();
//            real tDeltaT = tTimeCoeff.max() - tTimeCoeff.min();
//
//            // get tau2
//            real tTau2 = tDeltaT / 2;
//
//            // get tau3
//            real tTau3 = 0.25 * std::pow( tHugn, 2.0 ) / tPropConductivity->val()( 0 );
//
//            // compute stabilization parameter value
//            mPPVal = {{ std::pow(
//                    1 / std::pow( tTau1, 2.0 ) +
//                    1 / std::pow( tTau2, 2.0 ) +
//                    1 / std::pow( tTau3, 2.0 ) , -0.5 ) }};
        }

        //------------------------------------------------------------------------------
        void SP_SUPG_Advection::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the velocity FI
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // compute the diffusion coefficient
            std::shared_ptr< Property > tPropConductivity =
                    mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );
            real tK = tPropConductivity->val()( 0 );

            // evaluate norm( uTilde )
            real tNormA = std::sqrt( dot( tVelocityFI->val(), tVelocityFI->val() ) );

            // tau A
            real tTauA = 2.0 * tNormA / mElementSize;

            // tau K
            real tTauK = 4.0 * tK / std::pow( mElementSize, 2.0 );

            // init derivative of tauA, tauK, tauS
            Matrix< DDRMat > tdtauAdu( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );
            Matrix< DDRMat > tdtauKdu( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if dof type is velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity && tNormA > 0.0 )
            {
                // add contribution to dSPdu
                tdtauAdu.matrix_data() +=
                        2.0 * trans( tVelocityFI->val() ) * tVelocityFI->N() / ( mElementSize * tNormA );
            }

            // if conductivity property depends on dof type
            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                // compute tdtauKdu
                tdtauKdu.matrix_data() +=
                        4.0 * tPropConductivity->dPropdDOF( aDofTypes ) / std::pow( mElementSize, 2.0 );;
            }

            // evaluate tau
            Matrix< DDRMat > tTau = this->val();

            // scale dSPdu
            mdPPdMasterDof( tDofIndex ) = - 0.5 * std::pow( tTau( 0 ), 3 ) *
                    ( 2.0 * tTauA * tdtauAdu + 2.0 * tTauK * tdtauKdu );

//            // get the dof type index
//            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );
//
//            // get the dof type FI
//            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );
//
//            // set size for dSPdMasterDof, dTau1dDof, dTau3dDof
//            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//            Matrix< DDRMat > tdTau1dDof( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//            Matrix< DDRMat > tdTau3dDof( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//            Matrix< DDRMat > tdAlphadDof( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );
//
//            // get the velocity FI
//            Field_Interpolator * tVelocityFI
//            = mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );
//
//            // get the conductivity property
//            std::shared_ptr< Property > tPropConductivity
//            = mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );
//
//            // compute the norm of velocity
//            real tNorm = norm( tVelocityFI->val() );
//
//            // compute the abs term
//            real tAbs = 0.0;
//            uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
//            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
//            {
//                Matrix< DDRMat > tAdd = trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
//                tAbs += std::abs( tAdd( 0, 0 ) );
//            }
//
//            // compute hugn
//            real tHugn = 2.0 * tNorm / tAbs;
//
//            // compute tau1
//            real tTau1 = 0.5 * tHugn / tNorm;
//
//            // compute tau3
//            real tTau3 = 0.25 * std::pow( tHugn, 2.0 ) / tPropConductivity->val()( 0 );
//
//            // if dof type is velocity
//            if( aDofTypes( 0 ) == mMasterDofVelocity )
//            {
//                // compute derivative of the velocity norm
//                Matrix< DDRMat > tdNormdu = trans( tVelocityFI->val() ) * tVelocityFI->N() / tNorm;
//
//                // compute derivative of the abs term
//                Matrix< DDRMat > tdAbsdu( 1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );
//                uint tNumNodes = tVelocityFI->dnNdxn( 1 ).n_cols();
//                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
//                {
//                    Matrix< DDRMat > tAdd = trans( tVelocityFI->val() ) * tVelocityFI->dnNdxn( 1 ).get_column( iNode );
//                    tdAbsdu.matrix_data()
//                                    += tAdd * trans( tVelocityFI->dnNdxn( 1 ).get_column( iNode ) ) * tVelocityFI->N() / std::abs( tAdd( 0, 0 ) );
//                }
//
//                // compute derivative of hugn
//                Matrix< DDRMat > tdHugndu( 1, tVelocityFI->get_number_of_space_time_coefficients(), 0.0 );
//                tdHugndu = 2.0 * ( tdNormdu * tAbs - tdAbsdu * tNorm ) / std::pow( tAbs, 2.0 );
//
//                // compute dtau1du
//                tdTau1dDof.matrix_data() += 0.5 * ( tdHugndu * tNorm - tHugn * tdNormdu ) / std::pow( tNorm, 2.0 );
//
//                // compute dtau3du
//                tdTau3dDof.matrix_data() += 0.5 * tHugn * tdHugndu / tPropConductivity->val()( 0 );
//            }
//
//            // if conductivity property depends on dof type
//            if( tPropConductivity->check_dof_dependency( aDofTypes ) )
//            {
//                // compute dtau3du
//                tdTau3dDof.matrix_data() -=
//                        0.25 * std::pow( tHugn, 2.0 ) * tPropConductivity->dPropdDOF( aDofTypes ) / std::pow( tPropConductivity->val()( 0 ), 2.0 );
//            }
//
//            // add contribution
//            mdPPdMasterDof( tDofIndex ).matrix_data() +=
//                    std::pow( this->val()( 0 ), 3.0 ) *
//                    ( std::pow( tTau1, -3.0 ) * tdTau1dDof.matrix_data() + std::pow( tTau3, -3.0 ) * tdTau3dDof.matrix_data() );

        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


