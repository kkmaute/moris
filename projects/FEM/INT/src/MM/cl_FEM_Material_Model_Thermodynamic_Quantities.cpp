/*
 * cl_FEM_Material_Model_Thermodynamic_Quantities.cpp
 *
 *  Created on: Feb 12, 2021
 *      Author: wunsch
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

        const Matrix< DDRMat > & Material_Model::AlphaPDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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

        const Matrix< DDRMat > & Material_Model::BetaTDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
            // get the specific gas constant
            return get_property( "IsochoricHeatCapacity" )->val();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::CvDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
        const Matrix< DDRMat > & Material_Model::CvDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            return get_property( "IsochoricHeatCapacity" )->dPropdDOF( aDofTypes );
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
            // get the specific gas constant
            return get_property( "IsobaricHeatCapacity" )->val();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::CpDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
        const Matrix< DDRMat > & Material_Model::CpDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            return get_property( "IsobaricHeatCapacity" )->dPropdDOF( aDofTypes );
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
            // get the specific gas constant
            return get_property( "SpecificHeatRatio" )->val();
        }

        //------------------------------------------------------------------------------

        const Matrix< DDRMat > & Material_Model::GammaDOF_dep( const moris::Cell< MSI::Dof_Type > & aDofTypes )
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
        const Matrix< DDRMat > & Material_Model::GammaDOF_triv( const moris::Cell< MSI::Dof_Type > & aDofTypes )
        {
            // get the specific gas constant
            return get_property( "SpecificHeatRatio" )->dPropdDOF( aDofTypes );
        }             

        //-----------------------------------------------------------------------------

    }/* end_fem_namespace */
}/* end_moris_namespace */
