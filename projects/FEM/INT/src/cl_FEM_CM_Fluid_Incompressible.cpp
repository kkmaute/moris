
#include "cl_FEM_CM_Fluid_Incompressible.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "op_minus.hpp"

namespace moris
{
    namespace fem
    {

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::set_function_pointers()
        {
            switch ( mSpaceDim )
            {
                case ( 2 ):
                {
                    m_eval_strain       = &CM_Fluid_Incompressible::eval_strain_2d;
                    m_eval_divstrain    = &CM_Fluid_Incompressible::eval_divstrain_2d;
                    m_eval_teststrain   = &CM_Fluid_Incompressible::eval_teststrain_2d;
                    m_eval_dstraindx    = &CM_Fluid_Incompressible::eval_dstraindx_2d;
                    m_eval_ddivstraindu = &CM_Fluid_Incompressible::eval_ddivstraindu_2d;
                    m_flatten_normal    = &CM_Fluid_Incompressible::flatten_normal_2d;
                    break;
                }
                case ( 3 ):
                {
                    m_eval_strain       = &CM_Fluid_Incompressible::eval_strain_3d;
                    m_eval_divstrain    = &CM_Fluid_Incompressible::eval_divstrain_3d;
                    m_eval_teststrain   = &CM_Fluid_Incompressible::eval_teststrain_3d;
                    m_eval_dstraindx    = &CM_Fluid_Incompressible::eval_dstraindx_3d;
                    m_eval_ddivstraindu = &CM_Fluid_Incompressible::eval_ddivstraindu_3d;
                    m_flatten_normal    = &CM_Fluid_Incompressible::flatten_normal_3d;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "CM_Fluid_Incompressible::set_function_pointers - only works for 2d and 3d." );
                    break;
                }
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_flux()
        {
            // get the pressure FI
            Field_Interpolator * tPressureFI = mFIManager->get_field_interpolators_for_type( mDofMap[ "Pressure" ] );

            // create identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );

            // evaluate pressure contribution to flux
            Matrix< DDRMat > tP( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
            tP( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI * tPressureFI->val();

            // get the viscosity property
            std::shared_ptr< Property > tViscosityProp = mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // compute flux
            mFlux = -1.0 * tP + 2 * tViscosityProp->val()( 0 ) * this->strain();
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_divflux()
        {
            // get the pressure FI
            Field_Interpolator * tPressureFI = mFIManager->get_field_interpolators_for_type( mDofMap[ "Pressure" ] );

            // get the viscosity property
            std::shared_ptr< Property > tViscosityProp = mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // compute flux
            mDivFlux = -1.0 * tPressureFI->gradx( 1 ) + 2 * tViscosityProp->val()( 0 ) * this->divstrain();
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_ddivfluxdu( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the corresponding FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivflux/du
            mddivfluxdu( tDofIndex ).set_size( mSpaceDim, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity property
            std::shared_ptr< Property > tViscosityProp = mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofMap[ "Velocity" ] )
            {
                // fill ddivstrain/dv
                mddivfluxdu( tDofIndex ).matrix_data() += 2.0 * tViscosityProp->val()( 0 ) * this->ddivstraindu( aDofTypes );
            }

            // if pressure dof
            if( aDofTypes( 0 ) == mDofMap[ "Pressure" ] )
            {
                // fill ddivstrain/dp
                mddivfluxdu( tDofIndex ).matrix_data() -= tFI->dnNdxn( 1 ).matrix_data();
            }

            if( tViscosityProp->check_dof_dependency( aDofTypes ) )
            {
                // fill ddivstrain/du
                mddivfluxdu( tDofIndex ).matrix_data() += 2.0 * this->divstrain() * tViscosityProp->dPropdDOF( aDofTypes );
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dfluxdx( uint aOrder )
        {
            // only 1st order supported
            MORIS_ERROR( aOrder == 1, "CM_Fluid_Incompressible::eval_dfluxdx - only 1st order supported." );

            // get the pressure FI
            Field_Interpolator * tPressureFI
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Pressure" ] );

            // get the viscosity property
            std::shared_ptr< Property > tViscosityProp
            = mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // create identity matrix
            Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
            Matrix< DDRMat > tP( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
            tP( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

            // evaluate dfluxdx
            mdFluxdx( aOrder -1 ) = trans( tP ) * trans( tPressureFI->gradx( aOrder ) )
                                  - 2.0 * tViscosityProp->val() * this->dstraindx( aOrder );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_traction( const Matrix< DDRMat > & aNormal )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_traction - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_testTraction( const Matrix< DDRMat > & aNormal )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_testTraction - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_strain_2d()
        {
            // get the velocity spatial gradient from velocity FI
            Matrix< DDRMat > tVelocityGradx
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )
                        ->gradx( 1 );

            // evaluate the strain
            mStrain.set_size( 3, 1, 0.0 );
            mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
            mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
            mStrain( 2, 0 ) = tVelocityGradx( 1, 0 ) + tVelocityGradx( 0, 1 );
        }

        void CM_Fluid_Incompressible::eval_strain_3d()
        {
            // get the velocity spatial gradient from velocity FI
            Matrix< DDRMat > tVelocityGradx
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )
                        ->gradx( 1 );

            // evaluate the strain
            mStrain.set_size( 6, 1, 0.0 );
            mStrain( 0, 0 ) = tVelocityGradx( 0, 0 );
            mStrain( 1, 0 ) = tVelocityGradx( 1, 1 );
            mStrain( 2, 0 ) = tVelocityGradx( 2, 2 );
            mStrain( 3, 0 ) = tVelocityGradx( 1, 2 ) + tVelocityGradx( 2, 1 );
            mStrain( 4, 0 ) = tVelocityGradx( 0, 2 ) + tVelocityGradx( 2, 0 );
            mStrain( 5, 0 ) = tVelocityGradx( 0, 1 ) + tVelocityGradx( 1, 0 );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_divstrain_2d()
        {
            // set size for div strain
            mDivStrain.set_size( 2, 1, 0.0 );

            // get the velocity gradient
            Matrix< DDRMat > tVelocityGrad
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->gradx( 2 );

            // fill div strain
            mDivStrain( 0 ) = tVelocityGrad( 0, 0 )
                            + tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 );
            mDivStrain( 1 ) = tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 )
                            + tVelocityGrad( 1, 1 );
        }

        void CM_Fluid_Incompressible::eval_divstrain_3d()
        {
            // set size for div strain
            mDivStrain.set_size( 3, 1, 0.0 );

            // get the velocity gradient
            Matrix< DDRMat > tVelocityGrad
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->gradx( 2 );

            // fill div strain
            mDivStrain( 0 ) = tVelocityGrad( 0, 0 )
                            + tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 )
                            + tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 );
            mDivStrain( 1 ) = tVelocityGrad( 5, 0 ) + tVelocityGrad( 0, 1 )
                            + tVelocityGrad( 1, 1 )
                            + tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 );
            mDivStrain( 2 ) = tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 )
                            + tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 )
                            + tVelocityGrad( 2, 2 );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_ddivstraindu_2d( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 2, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if( aDofTypes( 0 ) == mDofMap[ "Velocity" ] )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // fill ddivstrain/du
                mddivstraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )             = tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 );
                mddivstraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) = tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )             = tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 );
            }
        }

        void CM_Fluid_Incompressible::eval_ddivstraindu_3d( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof type FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // set size for ddivstrain/du
            mddivstraindu( tDofIndex ).set_size( 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            if( aDofTypes( 0 ) == mDofMap[ "Velocity" ] )
            {
                // get the 2nd order derivative of the shape functions d2Ndx2
                Matrix< DDRMat > tVelocityd2Ndx2 = tFI->dnNdxn( 2 );

                // get number of space time bases
                uint tNumBases = tFI->get_number_of_space_time_bases();

                // fill ddivstrain/du
                mddivstraindu( tDofIndex )( { 0, 0 }, { 0, tNumBases - 1 } )                 = tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 ) + tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } )     = tVelocityd2Ndx2.get_row( 5 );
                mddivstraindu( tDofIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tVelocityd2Ndx2.get_row( 4 );

                mddivstraindu( tDofIndex )( { 1, 1 }, { 0, tNumBases - 1 } )                 = tVelocityd2Ndx2.get_row( 5 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 ) + tVelocityd2Ndx2.get_row( 2 );
                mddivstraindu( tDofIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tVelocityd2Ndx2.get_row( 3 );

                mddivstraindu( tDofIndex )( { 2, 2 }, { 0, tNumBases - 1 } )                 = tVelocityd2Ndx2.get_row( 4 );
                mddivstraindu( tDofIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } )     = tVelocityd2Ndx2.get_row( 3 );
                mddivstraindu( tDofIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tVelocityd2Ndx2.get_row( 0 ) + tVelocityd2Ndx2.get_row( 1 ) + tVelocityd2Ndx2.get_row( 2 );
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dstraindx_2d( uint aOrder )
        {
            switch ( aOrder )
            {
                case( 1 ) :
                {
                    // set size for dstraindx
                    mdStraindx( aOrder - 1 ).set_size( 3, 2, 0.0 );

                    // get the velocity gradient
                    Matrix< DDRMat > tVelocityGrad
                    = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->gradx( 2 );

                    // fill dstraindx
                    mdStraindx( aOrder - 1 )( 0, 0 ) = tVelocityGrad( 0, 0 );
                    mdStraindx( aOrder - 1 )( 0, 1 ) = tVelocityGrad( 2, 0 );
                    mdStraindx( aOrder - 1 )( 1, 0 ) = tVelocityGrad( 2, 1 );
                    mdStraindx( aOrder - 1 )( 1, 1 ) = tVelocityGrad( 1, 1 );
                    mdStraindx( aOrder - 1 )( 2, 0 ) = tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 );
                    mdStraindx( aOrder - 1 )( 2, 1 ) = tVelocityGrad( 1, 0 ) + tVelocityGrad( 2, 1 );
                }
                default :
                    MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dstraindx< 2 > - order not supported." );
            }
        }

        void CM_Fluid_Incompressible::eval_dstraindx_3d( uint aOrder )
        {
            switch ( aOrder )
            {
                case( 1 ) :
                {
                    // set size for dstraindx
                    mdStraindx( aOrder - 1 ).set_size( 6, 3, 0.0 );

                    // get the velocity gradient
                    Matrix< DDRMat > tVelocityGrad
                    = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->gradx( 2 );

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
                    mdStraindx( aOrder - 1 )( 3, 0 ) = tVelocityGrad( 4, 1 ) + tVelocityGrad( 5, 2 );
                    mdStraindx( aOrder - 1 )( 3, 1 ) = tVelocityGrad( 3, 1 ) + tVelocityGrad( 1, 2 );
                    mdStraindx( aOrder - 1 )( 3, 2 ) = tVelocityGrad( 2, 1 ) + tVelocityGrad( 3, 2 );
                    mdStraindx( aOrder - 1 )( 4, 0 ) = tVelocityGrad( 4, 0 ) + tVelocityGrad( 0, 2 );
                    mdStraindx( aOrder - 1 )( 4, 1 ) = tVelocityGrad( 3, 0 ) + tVelocityGrad( 5, 2 );
                    mdStraindx( aOrder - 1 )( 4, 2 ) = tVelocityGrad( 2, 0 ) + tVelocityGrad( 4, 2 );
                    mdStraindx( aOrder - 1 )( 5, 0 ) = tVelocityGrad( 2, 0 ) + tVelocityGrad( 0, 1 );
                    mdStraindx( aOrder - 1 )( 5, 1 ) = tVelocityGrad( 1, 0 ) + tVelocityGrad( 5, 1 );
                    mdStraindx( aOrder - 1 )( 5, 2 ) = tVelocityGrad( 3, 0 ) + tVelocityGrad( 4, 1 );
                }
                default :
                    MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dstraindx< 3 > - order not supported." );
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_teststrain_2d()
        {
            // compute displacement gradient
            Matrix< DDRMat > tdnNdxn
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 3, tNumBases * 2, 0.0 );
            mTestStrain( { 0, 0 },{ 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 },{ 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 },{ tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 2, 2 },{ tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

        void CM_Fluid_Incompressible::eval_teststrain_3d()
        {
            // compute displacement gradient
            Matrix< DDRMat > tdnNdxn
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->dnNdxn( 1 );

            // get number of bases for displacement
            uint tNumBases
            = mFIManager->get_field_interpolators_for_type( mDofMap[ "Velocity" ] )->get_number_of_space_time_bases();

            // build the test strain
            mTestStrain.set_size( 6, tNumBases * 3 , 0.0 );
            mTestStrain( { 0, 0 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
            mTestStrain( { 4, 4 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 5, 5 }, { 0, tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );

            mTestStrain( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );

            mTestStrain( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 2, 2 }, { 0, tNumBases - 1 } );
            mTestStrain( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 1, 1 }, { 0, tNumBases - 1 } );
            mTestStrain( { 4, 4 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = tdnNdxn( { 0, 0 }, { 0, tNumBases - 1 } );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dFluxdDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // get the dof FI
            Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init mdFluxdDof
            mdFluxdDof( tDofIndex ).set_size( ( mSpaceDim - 1) * 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // get the viscosity property
            std::shared_ptr< Property > tViscosityProp = mProperties( static_cast< uint >( CM_Property_Type::VISCOSITY ) );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofMap[ "Velocity" ] )
            {
                // build dfluxdv
                mdFluxdDof( tDofIndex ).matrix_data() += 2.0 * tViscosityProp->val()( 0 ) * this->dStraindDOF( aDofTypes );
            }

            // if pressure dof
            if ( aDofTypes( 0 ) == mDofMap[ "Pressure" ] )
            {
                // create identity matrix
                Matrix< DDRMat > tI( mSpaceDim, 1, 1.0 );
                Matrix< DDRMat > tII( ( mSpaceDim - 1 ) * 3, 1, 0.0 );
                tII( { 0, mSpaceDim - 1 }, { 0, 0 } ) = tI.matrix_data();

                // get shape function for presure field
                Matrix< DDRMat > tPressureN = tFI->N();

                // build the dfluxdp
                mdFluxdDof( tDofIndex ).matrix_data() -= tII * tPressureN;
            }

            // if indirect dependency on the dof type through viscosity
            if ( tViscosityProp->check_dof_dependency( aDofTypes ) )
            {
                // compute derivative with indirect dependency through properties
                mdFluxdDof( tDofIndex ).matrix_data() += 2.0 * this->strain() * tViscosityProp->dPropdDOF( aDofTypes );
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                          const Matrix< DDRMat >             & aNormal )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dTractiondDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dTestTractiondDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                              const Matrix< DDRMat >             & aNormal,
                                                              const Matrix< DDRMat >             & aJump )
        {
            MORIS_ERROR( false, "CM_Fluid_Incompressible::eval_dTestTractiondDOF - Not implemented." );
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::eval_dStraindDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the dof type as a uint
             uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

             // get the dof FI
             Field_Interpolator * tFI = mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

             // get the dof type index
             uint tDofIndex = mGlobalDofTypeMap( tDofType );

             // init mdStraindDof
             mdStraindDof( tDofIndex ).set_size( ( mSpaceDim - 1 ) * 3, tFI->get_number_of_space_time_coefficients(), 0.0 );

            // if velocity dof
            if( aDofTypes( 0 ) == mDofMap[ "Velocity" ] )
            {
                // compute derivative
                mdStraindDof( tDofIndex ).matrix_data() += this->testStrain().matrix_data();
            }
        }

//--------------------------------------------------------------------------------------------------------------
        void CM_Fluid_Incompressible::flatten_normal_2d( const Matrix< DDRMat > & aNormal,
                                                               Matrix< DDRMat > & aFlatNormal )
        {
            aFlatNormal.set_size( 2, 3, 0.0 );
            aFlatNormal( 0, 0 ) = aNormal( 0, 0 );
            aFlatNormal( 0, 2 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 1 ) = aNormal( 1, 0 );
            aFlatNormal( 1, 2 ) = aNormal( 0, 0 );
        }

        void CM_Fluid_Incompressible::flatten_normal_3d( const Matrix< DDRMat > & aNormal,
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

//--------------------------------------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
