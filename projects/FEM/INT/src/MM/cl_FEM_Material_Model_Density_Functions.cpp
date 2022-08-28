/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Material_Model_Density_Functions.cpp
 *
 */

#include "cl_FEM_Material_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        // RETURN FUNCTIONS FOR DENSITY (SECOND EQUATION OF STATE)
        //------------------------------------------------------------------------------

        // if thermodynamic state variable is dependent compute and retrieve value from storage
        const Matrix< DDRMat > & Material_Model::density_dep()
        {
            // if the density was not evaluated
            if( mDensityEval )
            {
                // evaluate the density
                this->eval_density();

                // set bool for evaluation
                mDensityEval = false;
            }
            // return the density value
            return mDensity;
        }

        // trivial operation: get value from FI
        const Matrix< DDRMat > & Material_Model::density_triv()
        {
            // return the density value
            return mFIManager->get_field_interpolators_for_type( mDofDensity )->val();
        }

        //------------------------------------------------------------------------------

        // if thermodynamic state variable is dependent compute and retrieve value from storage
        const Matrix< DDRMat > & Material_Model::DensityDot_dep()
        {
            // if the flux was not evaluated
            if( mDensityDotEval )
            {
                // evaluate the flux
                this->eval_DensityDot();

                // set bool for evaluation
                mDensityDotEval = false;
            }
            // return the flux value
            return mDensityDot;
        }

        // trivial operation: get value from FI
        const Matrix< DDRMat > & Material_Model::DensityDot_triv()
        {
            // return the density rate of change
            return mFIManager->get_field_interpolators_for_type( mDofDensity )->gradt( 1 );
        }

        //------------------------------------------------------------------------------

        // if thermodynamic state variable is dependent compute and retrieve value from storage
        const Matrix< DDRMat > & Material_Model::dnDensitydxn_dep( uint aOrder )
        {
            switch ( aOrder )
            {
                case 1: // first derivative
                {
                    if ( mdDensitydxEval )
                    {
                        // evaluate the flux
                        this->eval_dDensitydx();

                        // set bool for evaluation
                        mdDensitydxEval = false;
                    }

                    // return the flux value
                    return mdDensitydx;
                }

                case 2: // second derivative
                {

                    if ( md2Densitydx2Eval )
                    {
                        // evaluate the flux
                        this->eval_d2Densitydx2();

                        // set bool for evaluation
                        md2Densitydx2Eval = false;
                    }

                    // return the flux value
                    return md2Densitydx2;
                }

                default:
                {
                    MORIS_ERROR( false, "Material_Model::dnDensitydxn - aOrder unknown, only 1 and 2 supported." );
                    return mdDensitydx;
                }
            }
        }

        // trivial operation: get values from FI
        const Matrix< DDRMat > & Material_Model::dnDensitydxn_triv( uint aOrder )
        {
            // return the density rate of change
            return mFIManager->get_field_interpolators_for_type( mDofDensity )->gradx( aOrder );
        }

        //-----------------------------------------------------------------------------
        //-----------------------------------------------------------------------------

        // if thermodynamic state variable is dependent compute and retrieve values from storage
        const Matrix< DDRMat > & Material_Model::DensityDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the MM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofType ),
                    "Material_Model::DensityDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mDensityDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_DensityDOF( aDofType );

                // set bool for evaluation
                mDensityDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mDensityDof( tDofIndex );
        }

        // trivial operation: get values from FI
        const Matrix< DDRMat > & Material_Model::DensityDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // check DOF deriv is wrt to own DOF-type is with
            if ( aDofType( 0 ) != mDofDensity )
            {
                // get the dof type index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                if ( mDensityDofEval( tDofIndex ) )
                {
                    // initialize output matrix with zeros
                    mDensityDof( tDofIndex ).set_size( 1, // mSpaceDim
                            mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->
                            get_number_of_space_time_coefficients(), 0.0 );

                    // set flag
                    mDensityDofEval( tDofIndex )= false;
                }

                // return zero matrix
                return mDensityDof( tDofIndex );
            }
            else
            {
                // return the density rate of change
                return mFIManager->get_field_interpolators_for_type( mDofDensity )->N();
            }
        }

        //-----------------------------------------------------------------------------

        // if thermodynamic state variable is dependent compute and retrieve values from storage
        const Matrix< DDRMat > & Material_Model::DensityDotDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // if aDofType is not an active dof type for the MM
            MORIS_ASSERT(
                    this->check_dof_dependency( aDofType ),
                    "Material_Model::DensityDotDOF_dep - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mDensityDotDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_DensityDotDOF( aDofType );

                // set bool for evaluation
                mDensityDotDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mDensityDotDof( tDofIndex );
        }

        // trivial operation: get values from FI
        const Matrix< DDRMat > & Material_Model::DensityDotDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType )
        {
            // check DOF deriv is wrt to own DOF-type is with
            if ( aDofType( 0 ) != mDofDensity )
            {
                // get the dof type index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                if ( mDensityDotDofEval( tDofIndex ) )
                {
                    // initialize output matrix with zeros
                    mDensityDotDof( tDofIndex ).set_size( 1, // mSpaceDim
                            mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->
                            get_number_of_space_time_coefficients(), 0.0 );

                    // set flag
                    mDensityDotDofEval( tDofIndex )= false;
                }

                // return zero matrix
                return mDensityDotDof( tDofIndex );
            }
            else
            {
                // return the density rate of change
                return mFIManager->get_field_interpolators_for_type( mDofDensity )->dnNdtn( 1 );
            }
        }

        //-----------------------------------------------------------------------------

        // if thermodynamic state variable is dependent compute and retrieve values from storage
        const Matrix< DDRMat > & Material_Model::dnDensitydxnDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder )
        {
            // if aDofType is not an active dof type for the MM
            MORIS_ASSERT(
                    this->check_dof_dependency( aDofType ),
                    "Material_Model::dnDensitydxnDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

            switch ( aOrder )
            {
                case 1: // first derivative
                {
                    // if the derivative has not been evaluated yet
                    if( mdDensitydxDofEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_dDensitydxDOF( aDofType );

                        // set bool for evaluation
                        mdDensitydxDofEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return mdDensitydxDof( tDofIndex );
                }

                case 2: // second derivative
                {
                    // if the derivative has not been evaluated yet
                    if( md2Densitydx2DofEval( tDofIndex ) )
                    {
                        // evaluate the derivative
                        this->eval_d2Densitydx2DOF( aDofType );

                        // set bool for evaluation
                        md2Densitydx2DofEval( tDofIndex ) = false;
                    }

                    // return the derivative
                    return md2Densitydx2Dof( tDofIndex );
                }

                default:
                {
                    MORIS_ERROR( false, "Material_Model::dnDensitydxnDOF_dep - aOrder unknown, only 1 and 2 supported." );
                    return mdDensitydxDof( 0 );
                }
            }
        }

        // trivial operation: get values from FI
        const Matrix< DDRMat > & Material_Model::dnDensitydxnDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofType, uint aOrder )
        {
            // check DOF deriv is wrt to own DOF-type is with
            if ( aDofType( 0 ) != mDofDensity )
            {
                // get the dof type index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                if ( aOrder == 1 )
                {
                    if ( mdDensitydxDofEval( tDofIndex ) )
                    {
                        // initialize output matrix with zeros
                        mdDensitydxDof( tDofIndex ).set_size( mSpaceDim,
                                mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->
                                get_number_of_space_time_coefficients(), 0.0 );

                        // set flag
                        mdDensitydxDofEval( tDofIndex )= false;
                    }

                    // return zero matrix
                    return mdDensitydxDof( tDofIndex );
                }
                else if ( aOrder == 2 )
                {
                    if ( md2Densitydx2DofEval( tDofIndex ) )
                    {
                        // initialize output matrix with zeros
                        md2Densitydx2Dof( tDofIndex ).set_size( 3 * mSpaceDim - 3,
                                mFIManager->get_field_interpolators_for_type( aDofType( 0 ) )->
                                get_number_of_space_time_coefficients(), 0.0 );

                        // set flag
                        md2Densitydx2DofEval( tDofIndex )= false;
                    }

                    // return zero matrix
                    return md2Densitydx2Dof( tDofIndex );
                }
                else
                {
                    MORIS_ERROR( false, "Material_Model::dnDensitydxnDOF_triv - only orders 1 and 2 implemented." );
                    return mdDensitydxDof( 0 );
                }
            }
            else
            {
                // return the density rate of change
                return mFIManager->get_field_interpolators_for_type( mDofDensity )->dnNdxn( aOrder );
            }
        }

        //-----------------------------------------------------------------------------

    }/* end_fem_namespace */
}/* end_moris_namespace */

