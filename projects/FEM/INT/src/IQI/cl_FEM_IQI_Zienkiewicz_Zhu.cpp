/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Zienkiewicz_Zhu.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Zienkiewicz_Zhu.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Zienkiewicz_Zhu::IQI_Zienkiewicz_Zhu( enum Stress_Type aStressType )
    {
        // assign stress type to evaluate
        mStressType = aStressType;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Zienkiewicz_Zhu::compute_QI( Matrix< DDRMat >& aQI )
    {
        // initialize output
        aQI.set_size( 1, 1, 0.0 );

        // check if dof index was set
        if ( ( mLeaderDofTypes( 0 ).size() > 1 ) &&    //
                ( mStressType != Stress_Type::VON_MISES_STRESS ) )
        {
            MORIS_ERROR( mIQITypeIndex != -1,
                    "IQI_Zienkiewicz_Zhu::compute_QI - mIQITypeIndex not set, but is needed for principal, normal and shear stresses." );

            MORIS_ERROR( mIQITypeIndex <= 2 && mIQITypeIndex >= 0,
                    "IQI_Zienkiewicz_Zhu::compute_QI - mIQITypeIndex out of bounds, must be 0, 1, or 2 for the three spatial dimensions." );
        }

        // get field interpolator
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // initialize stress value
        real tStressValue = 0.0;

        // switch for different stress types
        switch ( mStressType )
        {
            case Stress_Type::VON_MISES_STRESS:
                tStressValue = this->eval_Von_Mises_stress();
                break;

            case Stress_Type::PRINCIPAL_STRESS:
                tStressValue = this->eval_principal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::NORMAL_STRESS:
                tStressValue = this->eval_normal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::SHEAR_STRESS:
                tStressValue = this->eval_shear_stress( mIQITypeIndex + 1 );
                break;

            default:
                MORIS_ERROR( false, "IQI_Zienkiewicz_Zhu::compute_QI - Unknown Stress Type." );
        }

        real tL2error = tFI->val()( 0 ) - tStressValue;

        // evaluate the QI
        aQI( 0, 0 ) = tL2error * tL2error;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Zienkiewicz_Zhu::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // check if dof index was set
        if ( ( mLeaderDofTypes( 0 ).size() > 1 ) &&    //
                ( mStressType != Stress_Type::VON_MISES_STRESS ) )
        {
            MORIS_ERROR( mIQITypeIndex != -1,
                    "IQI_Zienkiewicz_Zhu::compute_QI - mIQITypeIndex not set, but is needed for principal, normal and shear stresses." );

            MORIS_ERROR( mIQITypeIndex > 2,
                    "IQI_Zienkiewicz_Zhu::compute_QI - mIQITypeIndex out of bounds, must be 0, 1, or 2 for the three spatial dimensions." );
        }

        // get field interpolator
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // initialize stress value
        real tStressValue = 0.0;

        // switch for different stress types
        switch ( mStressType )
        {
            case Stress_Type::VON_MISES_STRESS:
                tStressValue = this->eval_Von_Mises_stress();
                break;

            case Stress_Type::PRINCIPAL_STRESS:
                tStressValue = this->eval_principal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::NORMAL_STRESS:
                tStressValue = this->eval_normal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::SHEAR_STRESS:
                tStressValue = this->eval_shear_stress( mIQITypeIndex + 1 );
                break;

            default:
                MORIS_ERROR( false, "IQI_Zienkiewicz_Zhu::compute_QI - Unknown Stress Type." );
        }

        real tL2error = tFI->val()( 0 ) - tStressValue;

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( tL2error * tL2error );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Zienkiewicz_Zhu::eval_Von_Mises_stress()
    {
        // get standardized stress vector
        Matrix< DDRMat > tStressVector;
        this->get_stress_vector( tStressVector );

        // compute contributions to von mises stress
        real tNormalStressContribution =
                std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 ) +    //
                std::pow( tStressVector( 1 ) - tStressVector( 2 ), 2.0 ) +    //
                std::pow( tStressVector( 2 ) - tStressVector( 0 ), 2.0 );
        real tShearStressContribution =
                std::pow( tStressVector( 3 ), 2.0 ) +    //
                std::pow( tStressVector( 4 ), 2.0 ) +    //
                std::pow( tStressVector( 5 ), 2.0 );

        real tVonMisesStress = std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );

        // compute Von-Mises stress value
        return tVonMisesStress;
    }

    //------------------------------------------------------------------------------

    real
    IQI_Zienkiewicz_Zhu::eval_principal_stress( uint aPrincipalStressIndex )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector;
        this->get_stress_vector( tStressVector );

        // initialize stress value
        real tStressValue = ( tStressVector( 0 ) + tStressVector( 1 ) ) / 2.0;

        // get max shear stress
        real tMaxShearStress = std::sqrt( std::pow( ( tStressVector( 0 ) - tStressVector( 1 ) ) / 2.0, 2.0 ) +    //
                                          std::pow( tStressVector( 2 ), 2.0 ) );

        // switch between the 2D case and the 3D case
        // for 2D
        if ( mLeaderFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions() == 2 )
        {
            switch ( aPrincipalStressIndex )
            {
                case 1:
                    tStressValue += tMaxShearStress;
                    break;

                case 2:
                    tStressValue -= tMaxShearStress;
                    break;

                default:
                    MORIS_ERROR( false,
                            "IQI_Zienkiewicz_Zhu::eval_principal_stress - Only 1st and 2nd principal stresses known for 2D." );
            }
        }

        // for 3D
        else
        {
            // compute invariants
            real tI1 = tStressVector( 0 ) + tStressVector( 1 ) + tStressVector( 2 );
            real tI2 = tStressVector( 0 ) * tStressVector( 1 ) +    //
                       tStressVector( 1 ) * tStressVector( 2 ) +    //
                       tStressVector( 2 ) * tStressVector( 0 ) -    //
                       std::pow( tStressVector( 3 ), 2.0 ) -        //
                       std::pow( tStressVector( 4 ), 2.0 ) -        //
                       std::pow( tStressVector( 5 ), 2.0 );
            real tI3 = tStressVector( 0 ) * tStressVector( 1 ) * tStressVector( 2 ) -    //
                       tStressVector( 0 ) * std::pow( tStressVector( 3 ), 2.0 ) -        //
                       tStressVector( 1 ) * std::pow( tStressVector( 4 ), 2.0 ) -        //
                       tStressVector( 3 ) * std::pow( tStressVector( 5 ), 2.0 ) +        //
                       2.0 * tStressVector( 3 ) * tStressVector( 4 ) * tStressVector( 5 );

            // help values
            real tQ = ( 1.0 / 9.0 ) * ( std::pow( tI1, 2.0 ) - 3.0 * tI2 );
            real tR = ( 1.0 / 54.0 ) * ( 2.0 * std::pow( tI1, 3.0 ) - 9.0 * tI1 * tI2 + 27.0 * tI3 );
            real tT = std::acos( tR / std::sqrt( std::pow( tQ, 3.0 ) ) );

            // compute principal values
            real tE1 = 2.0 * std::sqrt( tQ ) * std::cos( ( tT ) / 3.0 ) + tI1 / 3.0;
            real tE2 = 2.0 * std::sqrt( tQ ) * std::cos( ( tT + 2.0 * M_PI ) / 3.0 ) + tI1 / 3.0;
            real tE3 = 2.0 * std::sqrt( tQ ) * std::cos( ( tT + 4.0 * M_PI ) / 3.0 ) + tI1 / 3.0;

            // figure out minimum and maximum values
            switch ( aPrincipalStressIndex )
            {
                case 1:
                    tStressValue = std::max( std::max( tE1, tE2 ), std::max( tE2, tE3 ) );
                    break;

                case 2:
                    tStressValue = std::max( std::min( tE1, tE2 ), std::min( std::max( tE1, tE2 ), tE3 ) );
                    break;

                case 3:
                    tStressValue = std::min( std::min( tE1, tE2 ), std::min( tE2, tE3 ) );
                    break;

                default:
                    MORIS_ERROR( false,
                            "IQI_Zienkiewicz_Zhu::eval_principal_stress - Only 1st, 2nd, and 3rd principal stresses known for 3D." );
            }
        }

        // return stress value
        return tStressValue;
    }

    //------------------------------------------------------------------------------

    real
    IQI_Zienkiewicz_Zhu::eval_normal_stress( uint aStressIndex )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector;
        this->get_stress_vector( tStressVector );

        // pick stress value
        return tStressVector( aStressIndex - 1 );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Zienkiewicz_Zhu::eval_shear_stress( uint aStressIndex )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector;
        this->get_stress_vector( tStressVector );

        // pick stress value
        return tStressVector( aStressIndex + 2 );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Zienkiewicz_Zhu::get_stress_vector( Matrix< DDRMat >& aStressVector )
    {
        // create stress vector
        aStressVector.set_size( 6, 1, 0.0 );

        // get stress vector from Constitutive model
        Matrix< DDRMat > tCMStress =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux();

        // pull apart stress vector into components
        uint tNumStressComponents = tCMStress.length();
        switch ( tNumStressComponents )
        {
            // 2D plane stress
            case 3:
            {
                aStressVector( 0 ) = tCMStress( 0 );
                aStressVector( 1 ) = tCMStress( 1 );
                aStressVector( 5 ) = tCMStress( 2 );
                break;
            }
                // 2D plane strain and axisymmetric
            case 4:
            {
                aStressVector( 0 ) = tCMStress( 0 );
                aStressVector( 1 ) = tCMStress( 1 );
                aStressVector( 2 ) = tCMStress( 2 );
                aStressVector( 5 ) = tCMStress( 3 );
                break;
            }
                // 3D
            case 6:
            {
                aStressVector = tCMStress;
                break;
            }
                // Unknown size - error
            default:
                MORIS_ERROR( false,
                        "IQI_Zienkiewicz_Zhu::get_stress_vector - CM stress vector of unknown size; 3, 4 or 6 components expected." );
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
