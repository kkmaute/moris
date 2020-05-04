
//FEM/INT/src
#include "cl_FEM_SP_GGLS_Diffusion.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
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
            mPropertyMap[ "Conductivity" ] = Property_Type::CONDUCTIVITY;
            mPropertyMap[ "Density" ] = Property_Type::DENSITY;
            mPropertyMap[ "Heat_Capacity" ] = Property_Type::HEAT_CAPACITY;

            // populate the dof map (default)
            mMasterDofMap[ "Temp" ] = MSI::Dof_Type::TEMP;
        }

//------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::eval_SP()
        {
            // get the velocity FI
            Field_Interpolator * tVelocityFI
            = mMasterFIManager->get_field_interpolators_for_type( mMasterDofMap[ "Velocity" ] );

            // get the property values
            real tConductivity = mMasterProp( static_cast< uint >( Property_Type::CONDUCTIVITY ) )->val()( 0 );
            real tDensity      = mMasterProp( static_cast< uint >( Property_Type::DENSITY ) )->val()( 0 );
            real tHeatCapacity = mMasterProp( static_cast< uint >( Property_Type::HEAT_CAPACITY ) )->val()( 0 );

            // get alpha
            real tAlpha = (tDensity * tHeatCapacity / mTimStepSize) * std::pow(mElementSize, 2.0);

            // get xi-bar
            real tXiBar = ( std::cosh( std::sqrt(6*tAlpha) ) + 2 ) / ( std::cosh( std::sqrt(6*tAlpha) ) - 1 )  -  (1/tAlpha);

            // compute stabilization parameter value
            mPPVal = {{ ( std::pow(mElementSize, 2.0) / (6*tConductivity) ) * tXiBar }};
        }

//------------------------------------------------------------------------------
        void SP_GGLS_Diffusion::eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {

            MORIS_ERROR( false, "SP_SUPG_Advection::eval_dSPdMasterDOF - not implemented." );

        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


