
//FEM/INT/src
#include "cl_FEM_SP_GGLS_Diffusion.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_CM_Phase_State_Functions.hpp"
//LINALG/src
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        SP_GGLS_Diffusion::SP_GGLS_Diffusion()
        {
            // set the property pointer cell size
            mMasterProp.resize( static_cast< uint >( Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Conductivity" ]         = Property_Type::CONDUCTIVITY;
            mPropertyMap[ "Density" ]              = Property_Type::DENSITY;
            mPropertyMap[ "Heat_Capacity" ]        = Property_Type::HEAT_CAPACITY;
            mPropertyMap[ "Latent_Heat" ]          = Property_Type::LATENT_HEAT;
            mPropertyMap[ "PC_Temp" ]              = Property_Type::PC_TEMP;
            mPropertyMap[ "Phase_State_Function" ] = Property_Type::PHASE_STATE_FUNCTION;
            mPropertyMap[ "Phase_Change_Const" ]   = Property_Type::PHASE_CHANGE_CONST;
        }

        //------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::reset_cluster_measures()
        {
            // evaluate element size from the cluster
            mElementSize = mCluster->compute_cluster_cell_length_measure(
                    mtk::Primary_Void::PRIMARY,
                    mtk::Master_Slave::MASTER );
        }

        //------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                moris::Cell< std::string >                  & aDofStrings,
                mtk::Master_Slave                             aIsMaster )
        {
            // set dof type list
            mMasterDofTypes = aDofTypes;

            // switch on master slave
            switch ( aIsMaster )
            {
                case mtk::Master_Slave::MASTER :
                {
                    // loop on dof type
                    for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
                    {
                        // get dof string
                        std::string tDofString = aDofStrings( iDof );

                        // get dof type
                        MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                        // if velocity
                        if( "Temp" )
                        {
                            mMasterDofTemp = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    std::string( "SP_GGLS_Diffusion::set_dof_type_list - Unknown aDofString : ") +
                                    tDofString;
                            MORIS_ERROR( false , tErrMsg.c_str() );
                        }
                    }
                    break;
                }

                default:
                    MORIS_ERROR( false, "SP_GGLS_Diffusion::set_dof_type_list - unknown or incorrect master slave type." );
                    break;
            }
        }

        //------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR(
                    mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "SP_SUPG_Advection::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::eval_SP()
        {

            // get the property values
            real tConductivity = mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );
            real tDensity      = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            real tHeatCapacity = mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );
            std::shared_ptr< Property > tPropLatentHeat   = mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );


            // time step size
            real tDeltat = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();


            // get alpha
            real tAlpha = 0.0;

            // if there is a phase change
            if (tPropLatentHeat != nullptr)
            {
                // get values of properties
                real tLatentHeat = tPropLatentHeat->val()( 0 );
                real tMeltTemp   = mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) )->val()( 0 );
                real tPCconst    = mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) )->val()( 0 );
                real tPSfunc     = mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) )->val()( 0 );

                // get the dof type FI
                Field_Interpolator * tFIDer =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofTemp );

                // get phase state function
                moris::real tdfdT = eval_dFdTemp( tMeltTemp, tPCconst, tPSfunc, tFIDer);

                tAlpha = ( tDensity * (tHeatCapacity + tLatentHeat * tdfdT) * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }
            // if there is no phase change
            else
            {
                tAlpha = ( tDensity * tHeatCapacity * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }

            // compute stabilization parameter value
            mPPVal = {{ ( std::pow(mElementSize, 2.0) / 6.0 )
                    * ( ( std::cosh( std::sqrt(6.0*tAlpha) ) + 2 ) / ( std::cosh( std::sqrt(6.0*tAlpha) ) - 1 )  -  (1/tAlpha) )}};
        }


        //------------------------------------------------------------------------------
        // FIXME: this function needs separate testing
        void SP_GGLS_Diffusion::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {

            // get the property values
            std::shared_ptr< Property > tPropConductivity = mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) );
            std::shared_ptr< Property > tPropDensity      = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) );
            std::shared_ptr< Property > tPropHeatCapacity = mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) );

            std::shared_ptr< Property > tPropLatentHeat   = mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) );
            std::shared_ptr< Property > tPropMeltTemp     = mMasterProp( static_cast< uint >( Property_Type::PC_TEMP ) );
            std::shared_ptr< Property > tPropPCconst      = mMasterProp( static_cast< uint >( Property_Type::PHASE_CHANGE_CONST ) );
            std::shared_ptr< Property > tPropPSfunc       = mMasterProp( static_cast< uint >( Property_Type::PHASE_STATE_FUNCTION ) );


            // get the dof type index
            uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for dSPdMasterDof
            mdPPdMasterDof( tDofIndex ).set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );


            // time step size
            real tDeltat = mMasterFIManager->get_IP_geometry_interpolator()->get_time_step();

            // initialize values
            real tConductivity = tPropConductivity->val()( 0 );
            real tDensity      = tPropDensity->val()( 0 );
            real tHeatCapacity = tPropHeatCapacity->val()( 0 );
            real tLatentHeat = 0.0;
            real tMeltTemp   = 0.0;
            real tPCconst    = 0.0;
            real tPSfunc     = 0.0;
            real tAlpha      = 0.0;
            real tdfdT       = 0.0;

            // if there is a phase change
            if (tPropLatentHeat != nullptr)
            {

                // get values of properties
                tLatentHeat = tPropLatentHeat->val()( 0 );
                tMeltTemp   = tPropMeltTemp->val()( 0 );
                tPCconst    = tPropPCconst->val()( 0 );
                tPSfunc     = tPropPSfunc->val()( 0 );

                // get phase state function
                tdfdT = eval_dFdTemp( tMeltTemp, tPCconst, tPSfunc, tFIDer);

                tAlpha = ( tDensity * (tHeatCapacity + tLatentHeat * tdfdT) * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }

            // if there is no phase change
            else
            {
                tAlpha = ( tDensity * tHeatCapacity * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }

            // get derivative of SP wrt to alpha d(k*tau)/d(alpha)
            real dGammadAlpha =
                    -(   4 * std::cosh(std::sqrt(6*tAlpha))
            - 2 * std::pow(std::cosh(std::sqrt(6*tAlpha)), 2)
            + 3 * std::sqrt(6 * std::pow(tAlpha,3) ) * std::sinh(std::sqrt(6*tAlpha))
            - 2 ) /
            ( 2 * std::pow(tAlpha,2) * std::pow((std::cosh(std::sqrt(6*tAlpha)) - 1), 2) );

            // if dof type is termperature
            if( aDofTypes( 0 ) == mMasterDofTemp )
            {

                // if there is a phase change
                if (tPropLatentHeat != nullptr)
                {

                    // get Dof-deriv of phase state function
                    moris::Matrix<DDRMat> td2fdTdDOF = eval_dFdTempdDOF( tMeltTemp, tPCconst, tPSfunc, tFIDer);

                    mdPPdMasterDof( tDofIndex ).matrix_data() +=
                            dGammadAlpha * td2fdTdDOF * tDensity * tLatentHeat * std::pow(mElementSize, 2.0) / ( 6.0 * tConductivity * tDeltat );
                }
                // if no phase change defined, do nothing


            }

            // if there is a phase change
            if (tPropLatentHeat != nullptr)
            {
                // if indirect dependency of conductivity
                if ( mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() -=
                            dGammadAlpha * tPropConductivity->dPropdDOF( aDofTypes )
                            * (tHeatCapacity + tLatentHeat * tdfdT) * std::pow(mElementSize, 2.0)  / ( 6.0 * std::pow(tConductivity,2) * tDeltat );
                }

                // if indirect dependency of density
                if ( mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() +=
                            dGammadAlpha * tPropDensity->dPropdDOF( aDofTypes )
                            * (tHeatCapacity + tLatentHeat * tdfdT) * std::pow(mElementSize, 2.0)  / ( 6.0 * tConductivity * tDeltat );
                }

                // if indirect dependency of latent heat
                if ( mMasterProp( static_cast< uint >( Property_Type::LATENT_HEAT ) )->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() +=
                            dGammadAlpha * tDensity * tPropLatentHeat->dPropdDOF( aDofTypes ) * tdfdT
                            * std::pow(mElementSize, 2.0)  / ( 6.0 * tConductivity * tDeltat );
                }

            }

            // else if there is no phase change
            else
            {
                // if indirect dependency of conductivity
                if ( mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() -=
                            dGammadAlpha * tPropConductivity->dPropdDOF( aDofTypes )
                            * tHeatCapacity * std::pow(mElementSize, 2.0)  / ( 6.0 * std::pow(tConductivity,2) * tDeltat );
                }

                // if indirect dependency of density
                if ( mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() +=
                            dGammadAlpha * tPropDensity->dPropdDOF( aDofTypes )
                            * tHeatCapacity * std::pow(mElementSize, 2.0)  / ( 6.0 * tConductivity * tDeltat );
                }

            }

            // if indirect dependency of heat capacity
            if ( mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->check_dof_dependency( aDofTypes ) )
            {
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        dGammadAlpha * tDensity * tPropHeatCapacity->dPropdDOF( aDofTypes )
                        * std::pow(mElementSize, 2.0)  / ( 6.0 * tConductivity * tDeltat );
            }

        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


