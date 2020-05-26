
#include "cl_FEM_SP_Incompressible_Flow.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "fn_inv.hpp"
#include "op_div.hpp"
#include "fn_diag_vec.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------
        SP_Incompressible_Flow::SP_Incompressible_Flow()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Density" ]   = Property_Type::DENSITY;
            mPropertyMap[ "Viscosity" ] = Property_Type::VISCOSITY;
        }

        //--------------------------------------------------------------------------------------------------------------
        void SP_Incompressible_Flow::set_function_pointers()
        {
            // switch on space dimensions
            switch( mSpaceDim )
            {
                // if 2D
                case 2 :
                {
                    mEvalGFunc = this->eval_G_2d;
                    break;
                }

                // if 3D
                case 3 :
                {
                    mEvalGFunc = this->eval_G_3d;
                    break;
                }

                // default
                default :
                    MORIS_ERROR( false, "SP_Incompressible_Flow::set_function_pointers - only support 2 and 3D." );
                    break;
            }
        }

        //------------------------------------------------------------------------------
        void SP_Incompressible_Flow::set_dof_type_list(
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
                        if( "Velocity" )
                        {
                            mMasterDofVelocity = tDofType;
                        }
                        else if( "Pressure" )
                        {
                            mMasterDofPressure = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_Incompressible_Flow::set_dof_type_list - Unknown aDofString : ") +
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
                    MORIS_ERROR( false, "SP_Incompressible_Flow::set_dof_type_list - unknown master slave type." );
                    break;
            }
        }

        //------------------------------------------------------------------------------
        void SP_Incompressible_Flow::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_Incompressible_Flow::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_Incompressible_Flow::eval_SP()
        {
            // set size for SP values
            mPPVal.set_size( 2, 1, 0.0 );

            // get the velocity and pressure FIs
            Field_Interpolator * tVelocityFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the density and viscosity properties
            std::shared_ptr< Property > tDensityProp   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tViscosityProp = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // get the density and viscosity value
            real tDensity   = tDensityProp->val()( 0 );
            real tViscosity = tViscosityProp->val()( 0 );

            // get the parameter
            real tCI = mParameters( 0 )( 0 );

            // get the time step
            real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // evaluate Gij = sum_d dxi_d/dx_i dxi_d/dx_j
            Matrix< DDRMat > tG;
            this->eval_G( tG );

            // get flattened G
            Matrix< DDRMat > tFlatG = reshape( tG, 1, tG.numel() );

            // get trace of G
            real tTrG = sum( diag_vec( tG.matrix_data() ) );

            // evaluate tauM = mPPVal( 0 )
            Matrix< DDRMat > tvivjGij = trans( tVelocityFI->val() ) * tG * tVelocityFI->val();
            Matrix< DDRMat > tGijGij  = tFlatG * trans( tFlatG );
            real tPPVal = std::pow( 2.0 * tDensity / tDeltaT, 2.0 )
            + std::pow( tDensity, 2.0 ) * tvivjGij( 0 )
            + tCI * std::pow( tViscosity, 2.0 ) * tGijGij( 0 );
            mPPVal( 0 ) = std::pow( tPPVal, -0.5 );

            // evaluate tauC = mPPVal( 1 )
            mPPVal( 1 ) = 1.0 / ( mPPVal( 0 ) * tTrG );
        }

        //------------------------------------------------------------------------------
        void SP_Incompressible_Flow::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdPPdMasterDof( tDofIndex ).set_size( 2, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the velocity FI
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofVelocity );

            // get the density and viscosity properties
            std::shared_ptr< Property > tDensityProp   = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tViscosityProp = mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // get the density and viscosity value
            real tDensity   = tDensityProp->val()( 0 );
            real tViscosity = tViscosityProp->val()( 0 );

            // get the parameter
            real tCI = mParameters( 0 )( 0 );

            // get the time step
            real tDeltaT = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // evaluate Gij = sum_d dxi_d/dx_i dxi_d/dx_j
            Matrix< DDRMat > tG;
            this->eval_G( tG );

            // get flattened G
            Matrix< DDRMat > tFlatG = reshape( tG, 1, tG.numel() );

            // get trace of G
            real tTrG = sum( diag_vec( tG.matrix_data() ) );

            // evaluate
            Matrix< DDRMat > tvivjGij = trans( tVelocityFI->val() ) * tG * tVelocityFI->val();
            Matrix< DDRMat > tGijGij  = tFlatG * trans( tFlatG );

            // if velocity
            if( aDofTypes( 0 ) == mMasterDofVelocity )
            {
                mdPPdMasterDof( tDofIndex ).get_row( 0 ) +=
                        -0.5 * std::pow( this->val()( 0 ), 3 ) *
                        std::pow( tDensity, 2.0 ) *
                        ( 2.0 * trans( tFI->val() ) * tG * tFI->N() );
            }

            // if density
            if( tDensityProp->check_dof_dependency( aDofTypes ) )
            {
                mdPPdMasterDof( tDofIndex ).get_row( 0 ) +=
                        -0.5 * std::pow( this->val()( 0 ), 3 ) *
                        (   8.0 * tDensity / tDeltaT * tDensityProp->dPropdDOF( aDofTypes )
                                + 2.0 * tDensity * tvivjGij( 0 ) * tDensityProp->dPropdDOF( aDofTypes ) );
            }

            // if viscosity
            if( tViscosityProp->check_dof_dependency( aDofTypes ) )
            {
                mdPPdMasterDof( tDofIndex ).get_row( 0 ) +=
                        -0.5 * std::pow( this->val()( 0 ), 3 ) *
                        ( tCI * 2.0 * tViscosity * tGijGij * tViscosityProp->dPropdDOF( aDofTypes ) );
            }

            // dtauCdDOF
            mdPPdMasterDof( tDofIndex ).get_row( 1 ) = - 1.0 * mdPPdMasterDof( tDofIndex ).get_row( 0 ) / ( tTrG * std::pow( this->val()( 0 ), 2.0 ) );
        }

        //------------------------------------------------------------------------------
        void SP_Incompressible_Flow::eval_G( Matrix< DDRMat > & aG )
        {
            // get the space jacobian from IP geometry interpolator
            Matrix< DDRMat > tSpaceJacobian;
            mMasterFIManager->get_IP_geometry_interpolator()->space_jacobian( tSpaceJacobian );
            tSpaceJacobian = inv( tSpaceJacobian );

            // FIXME should not be here
            mSpaceDim = mMasterFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions();
            this->set_function_pointers();

            // evaluate Gij = sum_d dxi_d/dx_i dxi_d/dx_j
            this->mEvalGFunc( aG, tSpaceJacobian );
        }

        void SP_Incompressible_Flow::eval_G_2d( Matrix< DDRMat > & aG,
                Matrix< DDRMat > & aInvSpaceJacobian )
        {
            // set size for aG
            aG.set_size( 2, 2, 0.0 );

            // fill aGij = sum_d dxi_d/dx_i dxi_d/dx_j
            aG( 0, 0 ) = std::pow( aInvSpaceJacobian( 0, 0 ), 2.0 ) + std::pow( aInvSpaceJacobian( 0, 1 ), 2.0 );
            aG( 0, 1 ) = aInvSpaceJacobian( 0, 0 ) * aInvSpaceJacobian( 1, 0 ) + aInvSpaceJacobian( 0, 1 ) * aInvSpaceJacobian( 1, 1 );
            aG( 1, 0 ) = aG( 0, 1 );
            aG( 1, 1 ) = std::pow( aInvSpaceJacobian( 1, 0 ), 2.0 ) + std::pow( aInvSpaceJacobian( 1, 1 ), 2.0 );
        }

        void SP_Incompressible_Flow::eval_G_3d( Matrix< DDRMat > & aG,
                Matrix< DDRMat > & aInvSpaceJacobian )
        {
            // set size for aG
            aG.set_size( 3, 3, 0.0 );

            // fill aGij = sum_d dxi_d/dx_i dxi_d/dx_j
            aG( 0, 0 ) = std::pow( aInvSpaceJacobian( 0, 0 ), 2.0 )
            + std::pow( aInvSpaceJacobian( 0, 1 ), 2.0 )
            + std::pow( aInvSpaceJacobian( 0, 2 ), 2.0 );
            aG( 0, 1 ) = aInvSpaceJacobian( 0, 0 ) * aInvSpaceJacobian( 1, 0 )
                               + aInvSpaceJacobian( 0, 1 ) * aInvSpaceJacobian( 1, 1 )
                               + aInvSpaceJacobian( 0, 2 ) * aInvSpaceJacobian( 1, 2 );
            aG( 0, 2 ) = aInvSpaceJacobian( 0, 0 ) * aInvSpaceJacobian( 2, 0 )
                               + aInvSpaceJacobian( 0, 1 ) * aInvSpaceJacobian( 2, 1 )
                               + aInvSpaceJacobian( 0, 2 ) * aInvSpaceJacobian( 2, 2 );
            aG( 1, 0 ) = aG( 0, 1 );
            aG( 1, 1 ) = std::pow( aInvSpaceJacobian( 1, 0 ), 2.0 )
            + std::pow( aInvSpaceJacobian( 1, 1 ), 2.0 )
            + std::pow( aInvSpaceJacobian( 1, 2 ), 2.0 );
            aG( 1, 2 ) = aInvSpaceJacobian( 1, 0 ) * aInvSpaceJacobian( 2, 0 )
                               + aInvSpaceJacobian( 1, 1 ) * aInvSpaceJacobian( 2, 1 )
                               + aInvSpaceJacobian( 1, 2 ) * aInvSpaceJacobian( 2, 2 );
            aG( 2, 0 ) = aG( 0, 2 );
            aG( 2, 1 ) = aG( 1, 2 );
            aG( 2, 2 ) = std::pow( aInvSpaceJacobian( 2, 0 ), 2.0 )
            + std::pow( aInvSpaceJacobian( 2, 1 ), 2.0 )
            + std::pow( aInvSpaceJacobian( 2, 2 ), 2.0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


