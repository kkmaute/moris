//FEM/INT/src
#include "cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Spalart_Allmaras_Nitsche_Interface::SP_Spalart_Allmaras_Nitsche_Interface()
        {
            // set has slave flag to true
            mHasSlave = true;

            // set size for the property pointer cells
            mMasterProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mSlaveProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Material" ] = static_cast< uint >( SP_Property_Type::MATERIAL );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::reset_cluster_measures()
        {
            // evaluate cluster measures from the cluster
            mMasterVolume = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
            mSlaveVolume = mCluster->compute_cluster_cell_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::SLAVE );
            mInterfaceSurface = mCluster->compute_cluster_cell_side_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::set_dof_type_list(
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
                        if( tDofString == "Viscosity" )
                        {
                            mMasterDofViscosity = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_Spalart_Allmaras_Nitsche_Interface::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }

                case mtk::Master_Slave::SLAVE :
                {
                    // set dof type list
                    mSlaveDofTypes = aDofTypes;

                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( tDofString == "Viscosity" )
                        {
                            mSlaveDofViscosity = tDofType;
                        }
                        else
                        {
                            // error unknown dof string
                            MORIS_ERROR( false ,
                                    "SP_Spalart_Allmaras_Nitsche_Interface::set_dof_type_list - Unknown aDofString : %s \n",
                                    tDofString.c_str() );
                        }
                    }
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_Spalart_Allmaras_Nitsche_Interface::set_dof_type_list - unknown master slave type." );
            }
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_SP()
        {
            // init SP values
            mPPVal.set_size( 3, 1, 0.0 );

            // get the master/slave material diffusion coefficient values
            real tMasterMaterial = this->compute_diffusion_coefficient( mtk::Master_Slave::MASTER );
            real tSlaveMaterial  = this->compute_diffusion_coefficient( mtk::Master_Slave::SLAVE );

            // compute weighted property
            real tDeno = mMasterVolume / tMasterMaterial + mSlaveVolume / tSlaveMaterial;

            // compute Nitsche parameter value
            mPPVal( 0 ) = 2.0 * mParameters( 0 )( 0 ) * mInterfaceSurface / tDeno;

            // compute master weight value
            mPPVal( 1 ) = ( mMasterVolume / tMasterMaterial ) / tDeno;

            // compute slave weight value
            mPPVal( 2 ) = ( mSlaveVolume / tSlaveMaterial ) / tDeno;
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_dSPdMasterDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdMasterDof( tDofIndex ).set_size(
                    3,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the master/slave material diffusion coefficient values
            real tMasterMaterial = this->compute_diffusion_coefficient( mtk::Master_Slave::MASTER );
            real tSlaveMaterial  = this->compute_diffusion_coefficient( mtk::Master_Slave::SLAVE );

            // compute weighted property
            real tDeno = mMasterVolume / tMasterMaterial + mSlaveVolume / tSlaveMaterial;

            // get the master material diffusion coefficient derivatives
            Matrix< DDRMat > tMasterdmaterialdu;
            this->compute_ddiffusiondu( aDofTypes, tMasterdmaterialdu, mtk::Master_Slave::MASTER );

            // compute derivative of Nitsche SP FIXME ????
            mdPPdMasterDof( tDofIndex ).get_row( 0 ) +=
                    this->val()( 0 ) * mMasterVolume * tMasterdmaterialdu /
                    ( std::pow( tMasterMaterial, 2 ) * tDeno );

            // compute derivative of master weight SP
            mdPPdMasterDof( tDofIndex ).get_row( 1 ) -=
                    this->val()( 1 ) * mSlaveVolume * tMasterdmaterialdu /
                    ( tMasterMaterial * tSlaveMaterial * tDeno );

            // compute derivative of the slave weight SP
            mdPPdMasterDof( tDofIndex ).get_row( 2 ) +=
                    this->val()( 2 ) * mMasterVolume * tMasterdmaterialdu /
                    ( std::pow( tMasterMaterial, 2 ) * tDeno );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::eval_dSPdSlaveDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mSlaveGlobalDofTypeMap( tDofType );

            // get the derivative dof type FI
            Field_Interpolator * tFIDer =
                    mSlaveFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // reset the matrix
            mdPPdSlaveDof( tDofIndex ).set_size( 3, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the master/slave material diffusion coefficient values
            real tMasterMaterial = this->compute_diffusion_coefficient( mtk::Master_Slave::MASTER );
            real tSlaveMaterial  = this->compute_diffusion_coefficient( mtk::Master_Slave::SLAVE );

            // compute weighted property
            real tDeno = mMasterVolume / tMasterMaterial + mSlaveVolume / tSlaveMaterial;

            // get the slave material diffusion coefficient derivatives
            Matrix< DDRMat > tSlavedmaterialdu;
            this->compute_ddiffusiondu( aDofTypes, tSlavedmaterialdu, mtk::Master_Slave::SLAVE );

            // compute derivative of Nitsche SP FIXME ????
            mdPPdSlaveDof( tDofIndex ).get_row( 0 ) +=
                    this->val()( 0 ) * mSlaveVolume * tSlavedmaterialdu /
                    ( std::pow( tSlaveMaterial, 2 ) * tDeno );

            // compute derivative of master weight SP
            mdPPdSlaveDof( tDofIndex ).get_row( 1 ) +=
                    this->val()( 1 ) * mSlaveVolume * tSlavedmaterialdu /
                    ( std::pow( tSlaveMaterial, 2 ) * tDeno );

            // compute derivative of the slave weight SP
            mdPPdSlaveDof( tDofIndex ).get_row( 2 ) -=
                    this->val()( 2 ) * mMasterVolume * tSlavedmaterialdu /
                    ( tMasterMaterial * tSlaveMaterial * tDeno );
        }

        //------------------------------------------------------------------------------

        real SP_Spalart_Allmaras_Nitsche_Interface::compute_diffusion_coefficient(
                mtk::Master_Slave aIsMaster )
        {
            // init diffusion coeff
            real tDiffusionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mMasterDofViscosity );

            // get the wall distance property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the fluid kinematic viscosity value
            real tKinViscosity = tPropKinViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity ) / mSigma;
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = this->compute_fn( aIsMaster );

                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity * tFn ) / mSigma;
            }

            return tDiffusionTerm;
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::compute_ddiffusiondu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondu,
                mtk::Master_Slave                    aIsMaster )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( tDerDofType );

            // init ddiffusiondu
            addiffusiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mMasterDofViscosity );

            // get the kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // if derivative dof type is viscosity
                if( tDerDofType == mMasterDofViscosity )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = this->compute_fn( aIsMaster );

                // compute dfndu
                Matrix< DDRMat > tdfndu;
                this->compute_dfndu( aDofTypes, tdfndu, aIsMaster );

                // if derivative dof type is viscosity
                if( tDerDofType == mMasterDofViscosity )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tFn * tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }

                // add contribution from fn to ddiffusiondu
                addiffusiondu += tModViscosity * tdfndu / mSigma;
            }
        }

        //------------------------------------------------------------------------------

        real SP_Spalart_Allmaras_Nitsche_Interface::compute_fn(
                mtk::Master_Slave aIsMaster )
        {
            // compute chi, chi³
            real tChi  = this->compute_chi( aIsMaster );
            real tChi3 = std::pow( tChi, 3 );

            // compute fn
            return ( mCn1 + tChi3 ) / ( mCn1 - tChi3 );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::compute_dfndu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfndu,
                mtk::Master_Slave                    aIsMaster )
        {
            // compute chi, chi², chi³
            real tChi = this->compute_chi( aIsMaster );
            real tChi2 = std::pow( tChi, 2 );
            real tChi3 = std::pow( tChi, 3 );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu, aIsMaster );

            // compute adfndu
            adfndu = 6.0 * mCn1 * tChi2 * tdchidu / std::pow( mCn1 - tChi3, 2 );
        }

        //------------------------------------------------------------------------------

        real SP_Spalart_Allmaras_Nitsche_Interface::compute_chi(
                mtk::Master_Slave aIsMaster )
        {
            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mMasterDofViscosity );

            // get the kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // compute chi
            return tFIModViscosity->val()( 0 ) / tPropKinViscosity->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void SP_Spalart_Allmaras_Nitsche_Interface::compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidu,
                mtk::Master_Slave                    aIsMaster )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    mSet->get_field_interpolator_manager( aIsMaster )->
                    get_field_interpolators_for_type( mMasterDofViscosity );

            // get the kinematic viscosity property
            std::shared_ptr< Property > tPropKinViscosity =
                    this->get_properties( aIsMaster )( static_cast< uint >( SP_Property_Type::MATERIAL ) );

            // compute chi
            real tChi = tFIModViscosity->val()( 0 ) / tPropKinViscosity->val()( 0 );

            // if residual dof type (here viscosity)
            if( aDofTypes( 0 ) == mMasterDofViscosity )
            {
                adchidu += tDerFI->N() / tPropKinViscosity->val()( 0 );
            }

            if( tPropKinViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu -=
                        tChi * tPropKinViscosity->dPropdDOF( aDofTypes ) /
                        tPropKinViscosity->val()( 0 );
            }
        }


        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

