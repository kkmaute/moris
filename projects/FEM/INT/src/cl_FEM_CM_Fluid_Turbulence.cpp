//FEM/INt/src
#include "cl_FEM_CM_Fluid_Turbulence.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//LINALG/src
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {

        //--------------------------------------------------------------------------------------------------------------
        CM_Fluid_Turbulence::CM_Fluid_Turbulence()
        {
            // set the property pointer cell size
            mProperties.resize( static_cast< uint >( CM_Property_Type::MAX_ENUM ), nullptr );

            // populate the map
            mPropertyMap[ "Viscosity" ] = CM_Property_Type::VISCOSITY;
            mPropertyMap[ "Density" ]   = CM_Property_Type::DENSITY;
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::set_function_pointers()
        {
            switch ( mSpaceDim )
            {
                case ( 2 ):
                    {
                    m_eval_strain       = &CM_Fluid_Turbulence::eval_strain_2d;
                    m_eval_divstrain    = &CM_Fluid_Turbulence::eval_divstrain_2d;
                    m_eval_teststrain   = &CM_Fluid_Turbulence::eval_teststrain_2d;
                    m_eval_dstraindx    = &CM_Fluid_Turbulence::eval_dstraindx_2d;
                    m_eval_ddivstraindu = &CM_Fluid_Turbulence::eval_ddivstraindu_2d;
                    m_flatten_normal    = &CM_Fluid_Turbulence::flatten_normal_2d;
                    break;
                    }
                case ( 3 ):
                    {
                    m_eval_strain       = &CM_Fluid_Turbulence::eval_strain_3d;
                    m_eval_divstrain    = &CM_Fluid_Turbulence::eval_divstrain_3d;
                    m_eval_teststrain   = &CM_Fluid_Turbulence::eval_teststrain_3d;
                    m_eval_dstraindx    = &CM_Fluid_Turbulence::eval_dstraindx_3d;
                    m_eval_ddivstraindu = &CM_Fluid_Turbulence::eval_ddivstraindu_3d;
                    m_flatten_normal    = &CM_Fluid_Turbulence::flatten_normal_3d;
                    break;
                    }
                default :
                {
                    MORIS_ERROR( false, "CM_Fluid_Turbulence::set_function_pointers - only works for 2d and 3d." );
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::set_dof_type_list(
                moris::Cell< moris::Cell< MSI::Dof_Type > > aDofTypes,
                moris::Cell< std::string >                  aDofStrings )
        {
            // set dof type list
            Constitutive_Model::set_dof_type_list( aDofTypes );

            // loop over the provided dof type
            for( uint iDof = 0; iDof < aDofTypes.size(); iDof++ )
            {
                // get dof type string
                std::string tDofString = aDofStrings( iDof );

                // get dof type
                MSI::Dof_Type tDofType = aDofTypes( iDof )( 0 );

                // switch on dof type string
                if( tDofString == "Velocity" )
                {
                    mDofVelocity = tDofType;
                }
                else if( tDofString == "Viscosity" )
                {
                    mDofViscosity = tDofType;
                }
                else
                {
                    std::string tErrMsg =
                            std::string( "CM_Fluid_Turbulence::set_dof_type_list - Unknown aDofString : ") +
                            tDofString;
                    MORIS_ERROR( false , tErrMsg.c_str() );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::set_property(
                std::shared_ptr< fem::Property > aProperty,
                std::string                      aPropertyString )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "CM_Fluid_Turbulence::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // set the property in the property cell
            mProperties( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_flux()
        {
            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // compute flux
            mFlux = 2.0 * tViscosityT * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_divflux()
        {
            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // get the gradx for the turbulence viscosity
            Matrix< DDRMat > tdViscosityTdx;
            this->compute_dviscositytdx( tdViscosityTdx );

            // flatten dviscositytdx
            Matrix< DDRMat > tdViscosityTdxFlat;
            this->flatten_normal( tdViscosityTdx, tdViscosityTdxFlat );

            // compute flux
            mDivFlux = 2.0 * tViscosityT * this->divstrain()
                     + 2.0 * tdViscosityTdxFlat * this->strain();
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size(
                    mSpaceDim,
                    tFI->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // if velocity dof
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // add contribution to ddivstrain/dv
                mddivfluxdu( tDofIndex ).matrix_data() +=
                        2.0 * tViscosityT * this->ddivstraindu( aDofTypes );
            }

            // get the turbulence viscosity
            Matrix< DDRMat > tdviscositytdu;
            this->compute_dviscositytdu( aDofTypes, tdviscositytdu );

                // add contribution to ddivstrain/du
                mddivfluxdu( tDofIndex ).matrix_data() +=
                        2.0 * this->divstrain() * tdviscositytdu;

                // get the gradx for the turbulence viscosity
                Matrix< DDRMat > tdViscosityTdx;
                this->compute_dviscositytdx( tdViscosityTdx );

                // flatten dviscositytdx
                Matrix< DDRMat > tdViscosityTdxFlat;
                this->flatten_normal( tdViscosityTdx, tdViscosityTdxFlat );

                //
                mddivfluxdu( tDofIndex ).matrix_data() +=
                        2.0 * tdViscosityTdxFlat * this->dStraindDOF( aDofTypes );

                // get the derivative of gradx for the turbulence viscosity wrt dof
                Matrix< DDRMat > tdViscosityTdxdu;
                this->compute_dviscositytdxdu( aDofTypes, tdViscosityTdxdu );

                Matrix< DDRMat > tStrain = this->strain();
                Matrix< DDRMat > tStrainFull;
                switch ( mSpaceDim )
                {
                    case 2:
                    {
                        tStrainFull = {
                                { tStrain( 0 ), tStrain( 2 ) },
                                { tStrain( 2 ), tStrain( 1 ) } };
                        break;
                    }

                    case 3:
                    {
                        tStrainFull = {
                                { tStrain( 0 ), tStrain( 5 ), tStrain( 4 ) },
                                { tStrain( 5 ), tStrain( 1 ), tStrain( 3 ) },
                                { tStrain( 4 ), tStrain( 3 ), tStrain( 2 ) }};
                        break;
                    }

                    default:
                        MORIS_ERROR( false, "CM_Fluid_Turbulence::eval_ddivfluxdu - only 2 or 3D" );
                        break;
                }
                //
                mddivfluxdu( tDofIndex ).matrix_data() +=
                        2.0 * tStrainFull * tdViscosityTdxdu;
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_dfluxdx( uint aOrder )
        {
            // only 1st order supported
            MORIS_ERROR( aOrder == 1, "CM_Fluid_Incompressible::eval_dfluxdx - only 1st order supported." );

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // evaluate dfluxdx
            mdFluxdx( aOrder - 1 ) = - 2.0 * tViscosityT * this->dstraindx( aOrder );

            // FIXME need d/dx part related to tViscosityT
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute the traction
            mTraction = tFlatNormal * this->flux();
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_testTraction(
                const Matrix< DDRMat >             & aNormal,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the dof FI
            Field_Interpolator * tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // init mdFluxdDof
            mTestTraction( tTestDofIndex ).set_size(
                    mSpaceDim,
                    tFITest->get_number_of_space_time_coefficients(),
                    0.0 );

            // flatten the normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // if test traction wrt velocity
            if( aTestDofTypes( 0 ) == mDofVelocity )
            {
                // compute test traction wrt velocity
                mTestTraction( tTestDofIndex ).matrix_data() +=
                        2.0 * tViscosityT * tFlatNormal * this->testStrain();
            }

            // evaluate test turbulence viscosity
            Matrix< DDRMat > tdviscositytdutest;
            this->compute_dviscositytdu( aTestDofTypes, tdviscositytdutest );

            // add contribution from turbulence viscosity
            mTestTraction( tTestDofIndex ).matrix_data() +=
                    2.0 * tFlatNormal * this->strain() * tdviscositytdutest;

        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_strain_2d()
        {
            // get the velocity spatial gradient from velocity FI
            Matrix< DDRMat > tVelocityGradx =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

            // evaluate the strain
            mStrain.set_size( 3, 1, 0.0 );
            mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
            mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
            mStrain( 2, 0 ) = 0.5 * ( tVelocityGradx( 1, 0 ) + tVelocityGradx( 0, 1 ) );
        }

        void CM_Fluid_Turbulence::eval_strain_3d()
        {
            // get the velocity spatial gradient from velocity FI
            Matrix< DDRMat > tVelocityGradx =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 1 );

            // evaluate the strain
            mStrain.set_size( 6, 1, 0.0 );
            mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
            mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
            mStrain( 2, 0 ) = tVelocityGradx( 2, 2 );
            mStrain( 3, 0 ) = 0.5 * ( tVelocityGradx( 1, 2 ) + tVelocityGradx( 2, 1 ) );
            mStrain( 4, 0 ) = 0.5 * ( tVelocityGradx( 0, 2 ) + tVelocityGradx( 2, 0 ) );
            mStrain( 5, 0 ) = 0.5 * ( tVelocityGradx( 0, 1 ) + tVelocityGradx( 1, 0 ) );
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_divstrain_2d()
        {
            // set size for div strain
            mDivStrain.set_size( 2, 1, 0.0 );

            // get the velocity gradient
            Matrix< DDRMat > tVelocityGrad =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

            // fill div strain
            mDivStrain( 0 ) = tVelocityGrad( 0, 0 ) + 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 ) );
            mDivStrain( 1 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) ) + tVelocityGrad( 1, 1 );
        }

        void CM_Fluid_Turbulence::eval_divstrain_3d()
        {
            // set size for div strain
            mDivStrain.set_size( 3, 1, 0.0 );

            // get the velocity gradient
            Matrix< DDRMat > tVelocityGrad =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

            // fill div strain
            mDivStrain( 0 ) =
                    tVelocityGrad( 0, 0 ) +
                    0.5* ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 ) ) +
                    0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 ) );
            mDivStrain( 1 ) =
                    0.5 * ( tVelocityGrad( 5, 0 ) + tVelocityGrad( 0, 1 ) ) +
                    tVelocityGrad( 1, 1 ) +
                    0.5 * ( tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 ) );
            mDivStrain( 2 ) =
                    0.5 * ( tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 ) ) +
                    0.5 * ( tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 ) ) +
                    tVelocityGrad( 2, 2 );
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_ddivstraindu_2d( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size(
                    2,
                    tFI->get_number_of_space_time_coefficients(),
                    0.0 );

            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // fill ddivstrain/du
                mddivstraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )             = tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 );
                mddivstraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 );
            }
        }

        void CM_Fluid_Turbulence::eval_ddivstraindu_3d( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // fill ddivstrain/du
                mddivstraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )                 = tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 ) + 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 5 );
                mddivstraindu( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 4 );

                mddivstraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )                 = 0.5 * tVelocityd2Ndx2.get_row( 5 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 ) + 0.5 * tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 3 );

                mddivstraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )                 = 0.5 * tVelocityd2Ndx2.get_row( 4 );
                mddivstraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } )     = 0.5 * tVelocityd2Ndx2.get_row( 3 );
                mddivstraindu( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tVelocityd2Ndx2.get_row( 0 ) + 0.5 * tVelocityd2Ndx2.get_row( 1 ) + tVelocityd2Ndx2.get_row( 2 );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_dstraindx_2d( uint aOrder )
        {
            switch ( aOrder )
            {
                case( 1 ) :
                    {
                    // set size for dstraindx
                    mdStraindx( aOrder - 1 ).set_size( 3, 2, 0.0 );

                    // get the velocity gradient
                    Matrix< DDRMat > tVelocityGrad =
                            mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

                    // fill dstraindx
                    mdStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
                    mdStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 2, 0 );
                    mdStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 2, 1 );
                    mdStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
                    mdStraindx( aOrder - 1 )( 2, 0 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) );
                    mdStraindx( aOrder - 1 )( 2, 1 ) = 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 ) );
                    }
                default :
                    MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dstraindx_2d - order not supported." );
            }
        }

        void CM_Fluid_Turbulence::eval_dstraindx_3d( uint aOrder )
        {
            switch ( aOrder )
            {
                case( 1 ) :
                    {
                    // set size for dstraindx
                    mdStraindx( aOrder - 1 ).set_size( 6, 3, 0.0 );

                    // get the velocity gradient
                    Matrix< DDRMat > tVelocityGrad =
                            mFIManager->get_field_interpolators_for_type( mDofVelocity )->gradx( 2 );

                    // fill dstraindx
                    mdStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
                    mdStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 5, 0 );
                    mdStraindx( aOrder - 1 )( 0, 2 ) = tVelocityGrad( 4, 0 );
                    mdStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 5, 1 );
                    mdStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
                    mdStraindx( aOrder - 1 )( 1, 2 ) = tVelocityGrad( 3, 1 );
                    mdStraindx( aOrder - 1 )( 2, 0 ) = tVelocityGrad( 4, 2 );
                    mdStraindx( aOrder - 1 )( 2, 1 ) = tVelocityGrad( 3, 2 );
                    mdStraindx( aOrder - 1 )( 2, 2 ) = tVelocityGrad( 2, 2 );
                    mdStraindx( aOrder - 1 )( 3, 0 ) = 0.5 * ( tVelocityGrad( 4, 1 ) + tVelocityGrad( 5, 2 ) );
                    mdStraindx( aOrder - 1 )( 3, 1 ) = 0.5 * ( tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 ) );
                    mdStraindx( aOrder - 1 )( 3, 2 ) = 0.5 * ( tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 ) );
                    mdStraindx( aOrder - 1 )( 4, 0 ) = 0.5 * ( tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 ) );
                    mdStraindx( aOrder - 1 )( 4, 1 ) = 0.5 * ( tVelocityGrad( 3, 0 ) + tVelocityGrad( 5, 2 ) );
                    mdStraindx( aOrder - 1 )( 4, 2 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 ) );
                    mdStraindx( aOrder - 1 )( 5, 0 ) = 0.5 * ( tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 ) );
                    mdStraindx( aOrder - 1 )( 5, 1 ) = 0.5 * ( tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 ) );
                    mdStraindx( aOrder - 1 )( 5, 2 ) = 0.5 * ( tVelocityGrad( 3, 0 ) + tVelocityGrad( 4, 1 ) );
                    }
                default :
                    MORIS_ERROR( false, "CM_Fluid_Turbulence::eval_dstraindx_3d - order not supported." );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_teststrain_2d()
        {
            // compute displacement gradient
            Matrix< DDRMat > tdnNdxn =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 3, tNumBases * 2, 0.0 );
            mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        void CM_Fluid_Turbulence::eval_teststrain_3d()
        {
            // compute displacement gradient
            Matrix< DDRMat > tdnNdxn =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases =
                    mFIManager->get_field_interpolators_for_type( mDofVelocity )->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 6, tNumBases * 3, 0.0 );
            mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 4, 4 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 5, 5 }, { 0, tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

            mTestStrain( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = 0.5 * tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof FI
            Field_Interpolator * tFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            mdFluxdDof( tDofIndex ).set_size(
                    ( mSpaceDim - 1 ) * 3,
                    tFI->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the turbulence viscosity
            real tViscosityT = this->compute_turbulence_viscosity();

            // if velocity dof
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // build dfluxdv
                mdFluxdDof( tDofIndex ).matrix_data() +=
                        2.0 * tViscosityT * this->dStraindDOF( aDofTypes );
            }

            // evaluate derivative of the turbulence viscosity
            Matrix< DDRMat > tdviscositytdu;
            this->compute_dviscositytdu( aDofTypes, tdviscositytdu );

            // add contribution from turbulence viscosity
            mdFluxdDof( tDofIndex ).matrix_data() +=
                    2.0 * this->strain() * tdviscositytdu;
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_dTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute dtractiondu
            mdTractiondDof( tDofIndex ) = tFlatNormal * this->dFluxdDOF( aDofTypes );
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_dTestTractiondDOF(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                const Matrix< DDRMat >             & aNormal,
                const Matrix< DDRMat >             & aJump,
                const moris::Cell< MSI::Dof_Type > & aTestDofTypes )
        {
            // get test dof type index
            uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

            // get the test dof FI
            Field_Interpolator * tFITest =
                    mFIManager->get_field_interpolators_for_type( aTestDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init the dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).set_size(
                    tFITest->get_number_of_space_time_coefficients(),
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // evaluate derivative of the turbulence viscosity
            Matrix< DDRMat > tdviscositytduder;
            this->compute_dviscositytdu( aDofTypes, tdviscositytduder );

            // evaluate test of the turbulence viscosity
            Matrix< DDRMat > tdviscositytdutest;
            this->compute_dviscositytdu( aTestDofTypes, tdviscositytdutest );

            // flatten normal
            Matrix< DDRMat > tFlatNormal;
            this->flatten_normal( aNormal, tFlatNormal );

            // compute contribution to dTestTractiondDof
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ).matrix_data() +=
                    2.0 * trans( tFlatNormal * this->dStraindDOF( aTestDofTypes ) ) * aJump * tdviscositytduder +
                    2.0 * trans( tdviscositytdutest ) * trans( aJump ) * tFlatNormal * this->dStraindDOF( aDofTypes );

            // if viscosity is the test dof
            if( aTestDofTypes( 0 ) == mDofViscosity )
            {
                MORIS_ERROR( false, "CM_Fluid_Turbulence::eval_dTestTractiondDOF - Case not implemented so far, require d2viscositytdviscosity2" );
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
            uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

            // get the dof FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( tDofType );

            // init mdStraindDof
            mdStraindDof( tDofIndex ).set_size(
                    ( mSpaceDim - 1 ) * 3,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofVelocity )
            {
                // compute derivative
                mdStraindDof( tDofIndex ).matrix_data() += this->testStrain().matrix_data();
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::flatten_normal_2d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat >       & aFlatNormal )
        {
            aFlatNormal.set_size( 2, 3, 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
        }

        void CM_Fluid_Turbulence::flatten_normal_3d(
                const Matrix< DDRMat > & aNormal,
                Matrix< DDRMat > & aFlatNormal )
        {
            aFlatNormal.set_size( 3, 6, 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 2, 2 ) = aNormal( 2, 0 );
            aFlatNormal( 0, 4 ) = aNormal( 2, 0 );
            aFlatNormal( 0, 5 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 3 ) = aNormal( 2, 0 );
            aFlatNormal( 1, 5 ) = aNormal( 0, 0 );
            aFlatNormal( 2, 3 ) = aNormal( 1, 0 );
            aFlatNormal( 2, 4 ) = aNormal( 0, 0 );
        }


        //------------------------------------------------------------------------------
        real CM_Fluid_Turbulence::compute_turbulence_viscosity()
        {
            // init the turbulence viscosity
            real tViscosityT = 0.0;

            // get the viscosity dof type FI
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute turbulent viscosity
            tViscosityT = tFIViscosity->val()( 0 ) * tFv1;

            // return the turbulence viscosity
            return tViscosityT;
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dviscositytdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adviscositytdu )
        {
            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            adviscositytdu.set_size(
                    1,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the viscosity dof type FI
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mDofViscosity )
            {
                // compute fv1
                real tFv1 = this->compute_fv1();

                // add contribution to dSPdu
                adviscositytdu.matrix_data() += tFv1 * tFIViscosity->N();
            }

            // compute dfv1du
            Matrix< DDRMat > tdfv1du;
            this->compute_dfv1du( aDofTypes, tdfv1du );

            // add contribution from dfv1du
            adviscositytdu.matrix_data() += tFIViscosity->val() * tdfv1du;
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dviscositytdx( Matrix< DDRMat > & adviscositytdx )
        {
            // get the viscosity dof type FI
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // compute fv1
            real tFv1 = this->compute_fv1();

            // compute dfv1dx
            Matrix< DDRMat > tdfv1dx;
            this->compute_dfv1dx( tdfv1dx );

            // compute dviscositytdx
            adviscositytdx = tFIViscosity->gradx( 1 ) * tFv1 + tdfv1dx * tFIViscosity->val();
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dviscositytdxdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adviscositytdxdu )
        {
            // get the dof type FI
            Field_Interpolator * tFIDer =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set matrix size
            adviscositytdxdu.set_size(
                    mSpaceDim,
                    tFIDer->get_number_of_space_time_coefficients(),
                    0.0 );

            // get the viscosity dof type FI
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mDofViscosity )
            {
                // compute fv1
                real tFv1 = this->compute_fv1();

                // compute dfv1dx
                Matrix< DDRMat > tdfv1dx;
                this->compute_dfv1dx( tdfv1dx );

                // add contribution to dviscositytdxdu
                adviscositytdxdu.matrix_data() += tFv1 * tFIViscosity->dnNdxn( 1 ) + tdfv1dx * tFIViscosity->N();
            }

            // compute dfv1du
            Matrix< DDRMat > tdfv1du;
            this->compute_dfv1du( aDofTypes, tdfv1du );

            // compute dfv1dxdu
            Matrix< DDRMat > tdfv1dxdu;
            this->compute_dfv1dxdu( aDofTypes, tdfv1dxdu );

            // add contribution from dfv1du
            adviscositytdxdu.matrix_data() += tFIViscosity->gradx( 1 ) * tdfv1du + tFIViscosity->val()( 0 ) * tdfv1dxdu;
        }

        //------------------------------------------------------------------------------
        real CM_Fluid_Turbulence::compute_chi()
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // compute chi
            real tChi = tFIViscosity->val()( 0 ) / tPropViscosity->val()( 0 );

            return tChi;
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >             & adchidu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mDofViscosity )
            {
                adchidu.matrix_data() += tDerFI->N() / tPropViscosity->val()( 0 );
            }

            // if viscosity property depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidu.matrix_data() -=
                        tFIViscosity->val() * tPropViscosity->dPropdDOF( aDofTypes ) /
                        std::pow( tPropViscosity->val()( 0 ), 2 );
            }
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dchidx(Matrix< DDRMat > & adchidx )
        {
            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // compute dchidx
            adchidx = tFIViscosity->gradx( 1 ) / tPropViscosity->val()( 0 );

            // FIXME dependency of property on x not accounted for
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dchidxdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidxdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidxdu.set_size( mSpaceDim, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    mFIManager->get_field_interpolators_for_type( mDofViscosity );

            // get the density and gravity properties
            std::shared_ptr< Property > tPropViscosity =
                    mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == mDofViscosity )
            {
                adchidxdu.matrix_data() += tDerFI->dnNdxn( 1 ) / tPropViscosity->val()( 0 );
            }

            // if viscosity property depends on dof type
            if( tPropViscosity->check_dof_dependency( aDofTypes ) )
            {
                adchidxdu.matrix_data() -=
                        tFIViscosity->gradx( 1 ) * tPropViscosity->dPropdDOF( aDofTypes ) /
                        std::pow( tPropViscosity->val()( 0 ), 2 );
            }
        }

        //------------------------------------------------------------------------------
        real CM_Fluid_Turbulence::compute_fv1()
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute fv1
            real tFv1 = std::pow( tChi, 3.0 ) / ( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ) );

            return tFv1;
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dfv1du(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1du )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute adfv1du
            adfv1du =
                    3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu /
                    std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dfv1dx( Matrix< DDRMat > & adfv1dx )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidx
            Matrix< DDRMat > tdchidx;
            this->compute_dchidx( tdchidx );

            // compute dfv1dx
            adfv1dx =
                    3.0 * mCv1 * std::pow( tChi, 2 ) * tdchidx /
                    std::pow( std::pow( tChi, 3 ) + std::pow( mCv1, 3 ), 2 );

            // FIXME dependency of property on x not accounted for
        }

        //------------------------------------------------------------------------------
        void CM_Fluid_Turbulence::compute_dfv1dxdu(
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1dxdu )
        {
            // compute chi
            real tChi = this->compute_chi();

            // compute dchidx
            Matrix< DDRMat > tdchidx;
            this->compute_dchidx( tdchidx );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            this->compute_dchidu( aDofTypes, tdchidu );

            // compute dchidxdu
            Matrix< DDRMat > tdchidxdu;
            this->compute_dchidxdu( aDofTypes, tdchidxdu );

            // compute dfv1dxdu
            adfv1dxdu = 3.0 * mCv1 * (
                    2.0 * tChi * ( std::pow( mCv1, 3 ) - 2.0 * std::pow( tChi, 3 ) ) * tdchidx * tdchidu
                    + std::pow( tChi, 2 ) * ( std::pow( tChi, 3 ) + std::pow( mCv1, 3 ) ) * tdchidxdu ) /
                            std::pow( std::pow( tChi, 3 ) + std::pow( mCv1, 3 ), 3 );
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
