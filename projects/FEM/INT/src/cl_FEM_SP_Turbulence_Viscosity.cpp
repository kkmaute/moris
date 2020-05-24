//FEM/INT/src
#include "cl_FEM_SP_Turbulence_Viscosity.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------
        SP_Turbulence_Viscosity::SP_Turbulence_Viscosity()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ]    = Property_Type::VISCOSITY;

            // populate the dof map (default)
            mMasterDofMap[ "Viscosity" ] = MSI::Dof_Type::VISCOSITY;
        }

        //--------------------------------------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_Turbulence_Viscosity::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //--------------------------------------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::set_function_pointers()
        {
        }

        //------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::reset_cluster_measures()
        {
        }

        //------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::eval_SP()
        {
            // get the viscosity dof type FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Viscosity" ] );

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute turbulent viscosity
            mPPVal = tFIViscosity->val() * tFv1;
        }

        //------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity dof type FI
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Viscosity" ] );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mMasterDofMap[ "Viscosity" ] )
            {
                // compute fv1
                real tFv1 = this->compute_fv1();

                // add contribution to dSPdu
                mdPPdMasterDof( tDofIndex ).matrix_data() += tFv1 * tFIViscosity->N();
            }

            // compute dfv1du
            Matrix< DDRMat > tdfv1du;
            this->compute_dfv1du( aDofTypes, tdfv1du );

            // add contribution from dfv1du
            mdPPdMasterDof( tDofIndex ).matrix_data() += tFIViscosity->val() * tdfv1du;
        }

        //------------------------------------------------------------------------------
        real SP_Turbulence_Viscosity::compute_chi()
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Viscosity" ] );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );

            return tChi;
        }

        //------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::compute_dchidu(
                moris::Cell< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & adchidu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Viscosity" ] );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mMasterProp( static_cast< uint >( Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mMasterDofMap[ "Viscosity" ] )
            {
                adchidu.matrix_data() += tDerFI->N() / tPropViscosity->val()( 0 );
            }

            // if viscosity property depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu.matrix_data() -= tChi * tPropViscosity->dPropdDOF( aDofTypes ) / tPropViscosity->val()( 0 );
            }
        }

        //------------------------------------------------------------------------------
        real SP_Turbulence_Viscosity::compute_fv1()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute fv1
            real tFv1 = std::pow( tChi, 3.0 ) / ( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ) );

            return tFv1;
        }

        //------------------------------------------------------------------------------
        void SP_Turbulence_Viscosity::compute_dfv1du(
                moris::Cell< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & adfv1du )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute adfv1du
            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu / std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


