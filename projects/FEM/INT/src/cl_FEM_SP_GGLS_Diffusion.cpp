
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
                        if( tDofString == "Temp" )
                        {
                            mMasterDofTemp = tDofType;
                        }
                        else
                        {
                            // create error message
                            std::string tErrMsg =
                                    "SP_GGLS_Diffusion::set_dof_type_list - Unknown aDofString : " +
                                    tDofString;

                            // error
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
                    MORIS_ERROR( false, "SP_GGLS_Diffusion::set_dof_type_list - unknown master slave type." );
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
            if ( mPropertyMap.find( aPropertyString ) == mPropertyMap.end() )
            {
                // create error message
                std::string tErrMsg =
                        "SP_GGLS_Diffusion::set_property - Unknown aPropertyString: " +
                        aPropertyString;

                // error
                MORIS_ERROR( false , tErrMsg.c_str() );
            }

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
                Field_Interpolator * tFITemp =
                        mMasterFIManager->get_field_interpolators_for_type( mMasterDofTemp );

                // get phase state function
                moris::real tdfdT = eval_dFdTemp( tMeltTemp, tPCconst, tPSfunc, tFITemp );

                tAlpha = ( tDensity * (tHeatCapacity + tLatentHeat * tdfdT) * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }
            // if there is no phase change
            else
            {
                tAlpha = ( tDensity * tHeatCapacity * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }

            // xi-bar value
            moris::real tXibar = ( std::cosh( std::sqrt(6.0*tAlpha) ) + 2.0 ) / ( std::cosh( std::sqrt(6.0*tAlpha) ) - 1.0 )  -  (1.0/tAlpha);

            // compute stabilization parameter value
            mPPVal = {{  std::pow(mElementSize, 2.0) / 6.0  * tXibar }};

            /* Note:
             * the artificial GGLS conductivity is returned as the stabilization parameter,
             * which is equal to k*tau_GGLS as defined in Alberto Pizzolato's Thesis,
             * or equal to just tau as defined in Franca's 1988 paper on GGLS
             */
        }


        //------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the properties
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
            real tLatentHeat   = 0.0;
            real tMeltTemp     = 0.0;
            real tPCconst      = 0.0;
            real tPSfunc       = 0.0;
            real tAlpha        = 0.0;
            real tdfdT         = 0.0;

            // if there is a phase change
            if (tPropLatentHeat != nullptr)
            {
                // get values of properties
                tLatentHeat = tPropLatentHeat->val()( 0 );
                tMeltTemp   = tPropMeltTemp->val()( 0 );
                tPCconst    = tPropPCconst->val()( 0 );
                tPSfunc     = tPropPSfunc->val()( 0 );

                // get phase state function
                tdfdT = eval_dFdTemp( tMeltTemp, tPCconst, tPSfunc, tFIDer );

                tAlpha = ( tDensity * ( tHeatCapacity + tLatentHeat * tdfdT ) * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }

            // if there is no phase change
            else
            {
                tAlpha = ( tDensity * tHeatCapacity * std::pow(mElementSize, 2.0) ) / ( 6.0 * tConductivity * tDeltat );
            }

            // get derivative of SP wrt to alpha d(k*tau)/d(alpha)
            real dXibardAlpha =
                    (4.0 -
                    8.0 * std::cosh( std::sqrt(6.0*tAlpha) ) +
                    4.0 * std::pow( std::cosh( std::sqrt(6.0*tAlpha) ) , 2.0 ) -
                    std::pow( 6.0*tAlpha , 1.5 ) * std::sinh( std::sqrt(6.0*tAlpha) ) )
                                    / ( 4.0 * std::pow( tAlpha , 2.0 ) * std::pow( ( std::cosh( std::sqrt(6.0*tAlpha) ) - 1.0 ) , 2.0 ) );

            moris::real tdSPdAlpha = std::pow(mElementSize, 2.0) / 6.0 * dXibardAlpha;


            // indirect contributions for both with and without phase change ---------------------------

            // if indirect dependency on conductivity
            if ( tPropConductivity->check_dof_dependency( aDofTypes ) )
            {
                mdPPdMasterDof( tDofIndex ).matrix_data() -=
                        tdSPdAlpha *
                        ( tDensity * tHeatCapacity * std::pow(mElementSize, 2.0) /
                                ( 6.0 * std::pow(tConductivity, 2.0) * tDeltat ) ) *
                                tPropConductivity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency on density
            if ( tPropDensity->check_dof_dependency( aDofTypes ) )
            {
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        tdSPdAlpha *
                        ( tHeatCapacity * std::pow(mElementSize, 2.0) /
                                ( 6.0 * tConductivity * tDeltat ) ) *
                                tPropDensity->dPropdDOF( aDofTypes );
            }

            // if indirect dependency on heat capacity
            if ( tPropHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                mdPPdMasterDof( tDofIndex ).matrix_data() +=
                        tdSPdAlpha *
                        ( tDensity * std::pow(mElementSize, 2.0) /
                                ( 6.0 * tConductivity * tDeltat ) ) *
                                tPropHeatCapacity->dPropdDOF( aDofTypes );
            }

            // if there is a phase change --------------------------------------------------------------
            if (tPropLatentHeat != nullptr)
            {

                // if dof type is temperature
                if( aDofTypes( 0 ) == mMasterDofTemp )
                {
                    // get Dof-deriv of phase state function
                    moris::Matrix<DDRMat> td2fdTdDOF =
                            eval_dFdTempdDOF( tMeltTemp, tPCconst, tPSfunc, tFIDer);

                    // derivative of tau wrt temperature DOFs
                    mdPPdMasterDof( tDofIndex ).matrix_data() +=
                            tdSPdAlpha *
                            ( tDensity * tLatentHeat * std::pow(mElementSize, 2.0) /
                                    ( 6.0 * tConductivity * tDeltat ) ) *
                                    td2fdTdDOF;
                }

                // if indirect dependency on conductivity
                if ( tPropConductivity->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() -=
                            tdSPdAlpha *
                            ( tDensity * tLatentHeat * tdfdT * std::pow(mElementSize, 2.0) /
                                    ( 6.0 * std::pow(tConductivity, 2.0) * tDeltat ) ) *
                                    tPropConductivity->dPropdDOF( aDofTypes );
                }

                // if indirect dependency on density
                if ( tPropDensity->check_dof_dependency( aDofTypes ) )
                {
                    mdPPdMasterDof( tDofIndex ).matrix_data() +=
                            tdSPdAlpha *
                            ( tLatentHeat * tdfdT * std::pow(mElementSize, 2.0) /
                                    ( 6.0 * tConductivity * tDeltat ) ) *
                                    tPropDensity->dPropdDOF( aDofTypes );
                }

            }

        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


// BACKUP -----------------------------------------------
// -1.0 * (4.0 * std::cosh(std::sqrt(6.0*tAlpha)) -
//         2.0 * std::pow(std::cosh(std::sqrt(6.0*tAlpha)), 2.0) +
//         3.0 * std::sqrt(6.0 * std::pow(tAlpha,3.0) ) * std::sinh(std::sqrt(6.0*tAlpha)) - 2.0 ) /
//         ( 2.0 * std::pow(tAlpha,2.0) * std::pow((std::cosh(std::sqrt(6.0*tAlpha)) - 1.0), 2.0) );

