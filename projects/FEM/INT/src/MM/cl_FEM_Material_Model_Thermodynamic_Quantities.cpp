/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Material_Model_Thermodynamic_Quantities.cpp
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
        // RETURN FUNCTIONS FOR VOLUME EXPANSIVITY (ALPHA_P) AND ISOTHERMAL COMPRESSIBLILITY (BETA_T)
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::AlphaP()
        {
            // if the density was not evaluated
            if( mAlphaPEval )
            {
                // evaluate the density
                this->eval_VolumeExpansivity();

                // set bool for evaluation
                mAlphaPEval = false;
            }
            // return the density value
            return mAlphaP;
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::BetaT()
        {
            // if the density was not evaluated
            if( mBetaTEval )
            {
                // evaluate the density
                this->eval_IsothermalCompressibility();

                // set bool for evaluation
                mBetaTEval = false;
            }
            // return the density value
            return mBetaT;
        }

        //-----------------------------------------------------------------------------
        //-----------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::AlphaPDOF( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // if aDofTypes is not an active dof type for the MM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "Material_Model::AlphaPDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mAlphaPDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_VolumeExpansivityDOF( aDofTypes );

                // set bool for evaluation
                mAlphaPDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mAlphaPDof( tDofIndex );
        }

        //-----------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::BetaTDOF( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // if aDofTypes is not an active dof type for the MM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "Material_Model::BetaTDOF - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mBetaTDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_IsothermalCompressibilityDOF( aDofTypes );

                // set bool for evaluation
                mBetaTDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mBetaTDof( tDofIndex );
        }

        //------------------------------------------------------------------------------
        // RETURN FUNCTIONS FOR CV, CP, AND GAMMA
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::Cv_dep()
        {
            // if the density was not evaluated
            if( mCvEval )
            {
                // evaluate the density
                this->eval_Cv();

                // set bool for evaluation
                mCvEval = false;
            }
            // return the density value
            return mCv;
        }

        // trivial operation: get values from property
        const Matrix< DDRMat > & Material_Model::Cv_triv()
        {
            // get the spedific isochoric heat capacity
            return get_property( "IsochoricHeatCapacity" )->val();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::CvDOF_dep( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // if aDofTypes is not an active dof type for the MM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "Material_Model::CvDOF_dep - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mCvDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_CvDOF( aDofTypes );

                // set bool for evaluation
                mCvDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mCvDof( tDofIndex );
        }

        // trivial operation: get values from property
        const Matrix< DDRMat > & Material_Model::CvDOF_triv( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific isochorich heat capacity
            const std::shared_ptr< Property > tPropIsochoricHeatCapacity = get_property( "IsochoricHeatCapacity" );

            // check DoF dependency
            if ( tPropIsochoricHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                return tPropIsochoricHeatCapacity->dPropdDOF( aDofTypes );
            }
            else // no dependency of property on DoF type
            {
                // get the dof index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

                // initialize zero matrix
                // FIXME: this is slow
                mCvDof( tDofIndex ).set_size( 1,
                                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                                get_number_of_space_time_coefficients(), 0.0 );

                // return zero matrix
                return mCvDof( tDofIndex );
            }
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::Cp_dep()
        {
            // if the density was not evaluated
            if( mCpEval )
            {
                // evaluate the density
                this->eval_Cp();

                // set bool for evaluation
                mCpEval = false;
            }
            // return the density value
            return mCp;
        }

        // trivial operation: get values from property
        const Matrix< DDRMat > & Material_Model::Cp_triv()
        {
            // get the specific isobaric heat capacity
            return get_property( "IsobaricHeatCapacity" )->val();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::CpDOF_dep( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // if aDofTypes is not an active dof type for the MM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "Material_Model::CpDOF_dep - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mCpDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_CpDOF( aDofTypes );

                // set bool for evaluation
                mCpDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mCpDof( tDofIndex );
        }

        // trivial operation: get values from property
        const Matrix< DDRMat > & Material_Model::CpDOF_triv( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific isobaric heat capcity
            const std::shared_ptr< Property > tPropIsobaricHeatCapacity = get_property( "IsobaricHeatCapacity" );

            // check DoF dependency
            if ( tPropIsobaricHeatCapacity->check_dof_dependency( aDofTypes ) )
            {
                return tPropIsobaricHeatCapacity->dPropdDOF( aDofTypes );
            }
            else // no dependency of property on DoF type
            {
                // get the dof index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

                // initialize zero matrix
                // FIXME: this is slow
                mCpDof( tDofIndex ).set_size( 1,
                                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                                get_number_of_space_time_coefficients(), 0.0 );

                // return zero matrix
                return mCpDof( tDofIndex );
            }
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::Gamma_dep()
        {
            // if the density was not evaluated
            if( mGammaEval )
            {
                // evaluate the density
                this->eval_Gamma();

                // set bool for evaluation
                mGammaEval = false;
            }
            // return the density value
            return mGamma;
        }

        // trivial operation: get values from property
        const Matrix< DDRMat > & Material_Model::Gamma_triv()
        {
            // get the ratio of the specific heat capacities
            return get_property( "SpecificHeatRatio" )->val();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::GammaDOF_dep( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // if aDofTypes is not an active dof type for the MM
            MORIS_ERROR(
                    this->check_dof_dependency( aDofTypes ),
                    "Material_Model::GammaDOF_dep - no dependency in this dof type." );

            // get the dof index
            uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

            // if the derivative has not been evaluated yet
            if( mGammaDofEval( tDofIndex ) )
            {
                // evaluate the derivative
                this->eval_GammaDOF( aDofTypes );

                // set bool for evaluation
                mGammaDofEval( tDofIndex ) = false;
            }

            // return the derivative
            return mGammaDof( tDofIndex );
        }

        // trivial operation: get values from property
        const Matrix< DDRMat > & Material_Model::GammaDOF_triv( const Vector< MSI::Dof_Type > & aDofTypes )
        {
            // get the ratio of the specific heat capacities
            const std::shared_ptr< Property > tPropGamma = get_property( "SpecificHeatRatio" );

            // check DoF dependency
            if ( tPropGamma->check_dof_dependency( aDofTypes ) )
            {
                return tPropGamma->dPropdDOF( aDofTypes );
            }
            else // no dependency of property on DoF type
            {
                // get the dof index
                uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

                // initialize zero matrix
                // FIXME: this is slow
                mGammaDof( tDofIndex ).set_size( 1,
                                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) )->
                                get_number_of_space_time_coefficients(), 0.0 );

                // return zero matrix
                return mGammaDof( tDofIndex );
            }
        }

        //-----------------------------------------------------------------------------

    }/* end_fem_namespace */
}/* end_moris_namespace */

