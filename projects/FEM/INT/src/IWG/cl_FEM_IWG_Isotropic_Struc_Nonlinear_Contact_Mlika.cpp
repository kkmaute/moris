/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Mlika.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Nonlinear_Contact_Mlika.hpp"
#include "armadillo"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Vector.hpp"
#include "fn_dot_Arma.hpp"
#include "fn_inv.hpp"
#include "fn_isfinite.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
#include "fn_trans.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <memory>
#include <utility>

namespace moris::fem
{
    // TODO @ff: REMOVE
    [[maybe_unused]] std::string pmat( const Matrix< DDRMat >& aMatrix )
    {
        std::stringstream tStream;
        // set the precision of the stream
        tStream.precision( 10 );

        for ( uint i = 0; i < aMatrix.numel(); i++ )
        {
            tStream << aMatrix( i );
            if ( i < aMatrix.numel() - 1 )
            {
                tStream << ",";
            }
        }
        return tStream.str();
    }

    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::IWG_Isotropic_Struc_Nonlinear_Contact_Mlika( sint aBeta )
            : mBeta( aBeta )      // sign for symmetric/unsymmetric Nitsche
            , mTheta( -mBeta )    // TODO @ff check that the conversion is correct... check what mBeta stands for!
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
        mPropertyMap[ "Gap" ]       = static_cast< uint >( IWG_Property_Type::GAP );

        // set size for the constitutive model pointer cell
        // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
    }

    //------------------------------------------------------------------------------

    /**
     * \brief Implementation of this IWG is mainly based on
     * 'Mlika, Rabii. “Nitsche Method for Frictional Contact and Self-Contact: Mathematical and Numerical Study.” PhD Thesis, Universite de Lyon, 2018. https://theses.hal.science/tel-02067118.'
     * @attention The implementation is leader-oriented (integration happens on leader side) while Mlika (and the literature in general) integrates on the follower side. Thus, the coordinates X will be on the leader side and Y will be on the follower side (in the literature, it is the other way around).
     * \param aWStar
     */
    void
    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

        // get leader index for residual dof type, indices for assembly
        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        uint const tDofIndex      = mSet->get_dof_index_for_type( tDisplDofTypes( 0 ), mtk::Leader_Follower::LEADER );
        uint const tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
        uint const tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );
        auto const tResRange      = std::make_pair( tResStartIndex, tResStopIndex );

        // get field interpolator for the residual dof type
        Field_Interpolator*    tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

        // get follower field interpolator for the residual dof type
        Field_Interpolator*    tFollowerDofs     = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tFollowerGeometry = mFollowerFIManager->get_IG_geometry_interpolator();

        // get user defined constitutive model, stabilization parameter and thickness property
        const std::shared_ptr< Constitutive_Model >&      tConstitutiveModel      = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
        const std::shared_ptr< Stabilization_Parameter >& tStabilizationParameter = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
        const std::shared_ptr< Property >&                tThicknessProperty      = mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tThicknessProperty != nullptr ) ? tThicknessProperty->val()( 0 ) : 1;

        ///// MLIKA
        const real tNitscheParam = tStabilizationParameter->val()( 0 );    // stabilization parameter gamma

        // utility variables
        const uint             tDim      = mNormal.numel();
        const Matrix< DDRMat > tIdentity = eye( tDim, tDim );

        // displacement of leader (L) and follower (F)
        const Matrix< DDRMat > tUL = tLeaderDofs->val();      // ( 2 x 1 )
        const Matrix< DDRMat > tUF = tFollowerDofs->val();    // ( 2 x 1 )

        // gradient of the displacement field
        const Matrix< DDRMat > tGrad_UL = tLeaderDofs->gradx( 1 );      // ( 2 x 2 )
        const Matrix< DDRMat > tGrad_UF = tFollowerDofs->gradx( 1 );    // ( 2 x 2 )

        // test functions of leader (L) and follower (F)
        const Matrix< DDRMat > tNL = tLeaderDofs->N();      // ( 2 x 8 )
        const Matrix< DDRMat > tNF = tFollowerDofs->N();    // ( 2 x 8 )

        // gradient of the test functions
        const Matrix< DDRMat > tGrad_NL = tLeaderDofs->dnNdxn( 1 );      // (2 x 4)
        const Matrix< DDRMat > tGrad_NF = tFollowerDofs->dnNdxn( 1 );    // (2 x 4)

        // coordinates of leader (L) and follower (F) in the reference configuration
        const Matrix< DDRMat > tXL = trans( tLeaderGeometry->valx() );      // ( 2 x 1 )
        const Matrix< DDRMat > tXF = trans( tFollowerGeometry->valx() );    // ( 2 x 1 )

        Matrix< DDRMat > const txL = tLeaderGeometry->valx_current( tLeaderDofs );
        Matrix< DDRMat > const txF = tFollowerGeometry->valx_current( tFollowerDofs );

        // compute the deformation gradients at the leader and follower side
        const Matrix< DDRMat > tFL = tIdentity + tGrad_UL;    // (2 x 2)
        const Matrix< DDRMat > tFF = tIdentity + tGrad_UF;    // (2 x 2)

        const Matrix< DDRMat > tNormalLeader = tLeaderGeometry->get_normal_current( tLeaderDofs );

        const Matrix< DDRMat > tNormalProjector = tNormalLeader * trans( tNormalLeader );

        // gap between points in current configuration
        const real             tInitialGap = dot( ( tXF - tXL ), tNormalLeader );
        const Matrix< DDRMat > tJump       = tUL - tUF - tNormalLeader * tInitialGap;

        // evaluate traction and their pressure equivalent (normal component)
        const Matrix< DDRMat > tTraction       = tConstitutiveModel->traction( tNormalLeader );    // ( 2 x 1 )
        const real             tPressure       = dot( tTraction, tNormalLeader );
        const Matrix< DDRMat > tNormalTraction = tPressure * tNormalLeader;    // ( 2 x 1 ) normal component of traction

        // evaluate test traction and the corresponding pressure equivalent (normal component) by computing the dot product between
        // each column of tTestTraction and tNormalLeader. This is done by element-wise multiplication and then summing the columns.
        const Matrix< DDRMat > tTestTraction       = tConstitutiveModel->testTraction( tNormalLeader, tDisplDofTypes );    // ( 2 x 8 )
        const Matrix< DDRMat > tTestPressure       = trans( tNormalLeader ) * tTestTraction;                               // ( 1 x 8 )
        const Matrix< DDRMat > tNormalTestTraction = tNormalLeader * tTestPressure;                                        // ( 2 x 8 ) normal component of test traction

        // Compute the contact term C_gamma(sigma, g, n), see Mlika (2018) Eq. (4.2.1). If the contact term is negative, the second term of the residual is not zero.
        const real tContactTerm = tPressure - tNitscheParam * dot( tJump, tNormalLeader );
        if ( tContactTerm < 0 )
        {
            mSet->get_residual()( 0 )( tResRange ) += aWStar * 0.5 * (                                                        //
                                                              -trans( tNL ) * tNormalProjector * tTraction                    //
                                                              - mTheta * trans( tTestTraction ) * tNormalProjector * tJump    //
                                                              + tNitscheParam * trans( tNL ) * tNormalProjector * tJump       //
                                                      );
        }
        else
        {
            // first integral in Mlika (2018) Eq. (4.25) (ensures 0 pressure when not in contact)
            mSet->get_residual()( 0 )( tResRange ) += aWStar * 0.5 * ( ( -mTheta / tNitscheParam ) * trans( tTestTraction ) * tNormalProjector * tTraction );
        }

        [[maybe_unused]] uint const tIteration = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve" );
        const Matrix< DDRMat >      tDebugArea{
                 { 0.55, 0.60 },    // xmin, xmax
                 { 0.5, 0.8 }       // ymin, ymax
        };
        // check if point tXL is in the debug area
        [[maybe_unused]] const bool tIsDebugPoint = ( tXL( 0 ) > tDebugArea( 0, 0 )
                                                      && tXL( 0 ) < tDebugArea( 0, 1 )
                                                      && tXL( 1 ) > tDebugArea( 1, 0 )
                                                      && tXL( 1 ) < tDebugArea( 1, 1 ) );
        if ( tIsDebugPoint )
        {
            std::cout << "IWG:" << std::setprecision( 10 )
                      << tIteration << ","
                      << dot( tJump, tNormalLeader ) << ","
                      << pmat( tXL ) << ","
                      << pmat( tXF ) << ","
                      << pmat( tUL ) << ","
                      << pmat( tUF ) << ","
                      << pmat( txL ) << ","
                      << pmat( txF ) << ","
                      << pmat( tNormalLeader ) << ","
                      << std::endl;
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                "IWG_Isotropic_Struc_Linear_Contact_Gap::compute_residual - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
        // check leader and follower field interpolators
        this->check_field_interpolators( mtk::Leader_Follower::LEADER );
        this->check_field_interpolators( mtk::Leader_Follower::FOLLOWER );
#endif

        Vector< MSI::Dof_Type > const tDisplDofTypes = mResidualDofType( 0 );

        // get leader index for residual dof type, indices for assembly
        uint const tDofIndex      = mSet->get_dof_index_for_type( tDisplDofTypes( 0 ), mtk::Leader_Follower::LEADER );
        uint const tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
        uint const tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );
        auto const tResRange      = std::make_pair( tResStartIndex, tResStopIndex );

        // get field interpolator for the residual dof type
        Field_Interpolator*    tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

        // get follower field interpolator for the residual dof type
        Field_Interpolator*    tFollowerDofs     = mFollowerFIManager->get_field_interpolators_for_type( tDisplDofTypes( 0 ) );
        Geometry_Interpolator* tFollowerGeometry = mFollowerFIManager->get_IG_geometry_interpolator();

        // get user defined constitutive model, stabilization parameter and thickness property
        std::shared_ptr< Constitutive_Model >&            tConstitutiveModel      = mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );
        const std::shared_ptr< Stabilization_Parameter >& tStabilizationParameter = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
        const std::shared_ptr< Property >&                tThicknessProperty      = mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tThicknessProperty != nullptr ) ? tThicknessProperty->val()( 0 ) : 1;

        const real tNitscheParam = tStabilizationParameter->val()( 0 );    // stabilization parameter gamma

        uint const tDim = mNormal.numel();

        // displacement of leader (L) and follower (F)
        const Matrix< DDRMat > tUL = tLeaderDofs->val();      // ( 2 x 1 )
        const Matrix< DDRMat > tUF = tFollowerDofs->val();    // ( 2 x 1 )

        // gradient of the displacement field
        const Matrix< DDRMat > tGrad_UL = tLeaderDofs->gradx( 1 );      // ( 2 x 2 )
        const Matrix< DDRMat > tGrad_UF = tFollowerDofs->gradx( 1 );    // ( 2 x 2 )

        // test functions of leader (L) and follower (F)
        // TODO @ff: Is this the correct representation of delta u? Can I be sure that the hat(delta u) always cancel out in the integrals?
        // TODO @ff: Why are the shape functions for x and y direction not evaluated in the same column? I.e. why the block structure with 0s?
        const Matrix< DDRMat > tdUL = tLeaderDofs->N();      // ( 2 x 8 )
        const Matrix< DDRMat > tdUF = tFollowerDofs->N();    // ( 2 x 8 )

        // gradient of the test functions
        const Matrix< DDRMat > tGrad_dUL = tLeaderDofs->dnNdxn( 1 );      // (2 x 4)
        const Matrix< DDRMat > tGrad_dUF = tFollowerDofs->dnNdxn( 1 );    // (2 x 4)

        // coordinates of leader (L) and follower (F) in the reference configuration
        const Matrix< DDRMat > tXL = trans( tLeaderGeometry->valx() );      // ( 2 x 1 )
        const Matrix< DDRMat > tXF = trans( tFollowerGeometry->valx() );    // ( 2 x 1 )

        // coordinates of leader (L) and follower (F) in the current configuration
        const Matrix< DDRMat > txL = tXL + tUL;    // ( 2 x 1 )
        const Matrix< DDRMat > txF = tXF + tUF;    // ( 2 x 1 )

        const Matrix< DDRMat > tIdentity = eye( tDim, tDim );

        // compute the deformation gradients at the leader and follower side
        Matrix< DDRMat > const tFL = tIdentity + tGrad_UL;    // (2 x 2)
        Matrix< DDRMat > const tFF = tIdentity + tGrad_UF;    // (2 x 2)

        // normal of leader (L) and follower (F) in the deformed configuration
        const Matrix< DDRMat > tNormalL = tFL * mNormal;
        const Matrix< DDRMat > tNormalF = tFF * tFollowerGeometry->get_normal();

        const Matrix< DDRMat > tNormalProjector = tNormalL * trans( tNormalL );

        // gap between points in current configuration
        const real tGap = dot( ( txF - txL ), tNormalL );

        // TODO: This Jump is not correct (dummy value)
        const real tJump = tGap;

        // evaluate traction and their pressure equivalent (normal component)
        const Matrix< DDRMat > tTraction       = tConstitutiveModel->traction( tNormalL, CM_Function_Type::PK1 );    // ( 2 x 1 )
        const real             tPressure       = dot( tTraction, tNormalL );
        const Matrix< DDRMat > tNormalTraction = tPressure * tNormalL;    // ( 2 x 1 ) normal component of traction

        // evaluate test traction and the corresponding pressure equivalent (normal component)
        const Matrix< DDRMat > tTestTraction       = tConstitutiveModel->testTraction( tNormalL, tDisplDofTypes, CM_Function_Type::PK1 );    // ( 2 x 8 )
        const Matrix< DDRMat > tTestPressure       = trans( tNormalL ) * tTestTraction;                                                      // ( 1 x 8 )
        const Matrix< DDRMat > tNormalTestTraction = tNormalL * tTestPressure;                                                               // ( 2 x 8 ) normal component of test traction

        // Compute the contact term C_gamma(sigma, g, n), see Mlika (2018) Eq. (4.2.1).
        const real tContactTerm = tPressure + tNitscheParam * tGap;

        // compute the Jacobian for all dof dependencies on the leader
        for ( uint iDOF = 0; iDOF < mRequestedLeaderGlobalDofTypes.size(); iDOF++ )
        {
            // get the dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

            // get the index for the dof type
            const sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            const uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
            const uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );
            const auto tDepRange      = std::make_pair( tDepStartIndex, tDepStopIndex );

            // extract sub-matrices
            auto& tJac = mSet->get_jacobian();

            // if constitutive model has a dependency on the dof type
            // if ( tConstitutiveModel->check_dof_dependency( tDofType ) )
            if ( tDofType( 0 ) == tDisplDofTypes( 0 ) )
            {
                const Matrix< DDRMat > tTraction_dU     = tConstitutiveModel->dTractiondDOF( tDofType, tNormalL, CM_Function_Type::PK1 );                                                  // ( 2 x 8 )
                const Matrix< DDRMat > tTestTraction_dU = tConstitutiveModel->dTestTractiondDOF( tDofType, tNormalL, tNormalProjector * tJump, tDisplDofTypes, CM_Function_Type::PK1 );    // ( 8 x 8 )

                // Compute derivatives of the contact term w.r.t the traction, the gap and the normal
                // All terms are zero by default and will be computed if the heaviside function of
                //    H(-tPressure - tNitscheParam * tGap) > 0
                // see Mlika (2018) Appendix B
                Matrix< DDRMat > tContactTerm_dTraction( tDim, tDim, 0.0 );    // ( 2 x 2 )
                Matrix< DDRMat > tContactTerm_dGap( tDim, 1, 0.0 );            // ( 2 x 1 )
                Matrix< DDRMat > tContactTerm_dNormal( tDim, tDim, 0.0 );      // ( 2 x 2 )
                if ( true || -tPressure - tNitscheParam * tGap > 0.0 )         // TODO @ff: remove true (for debugging)
                {
                    tContactTerm_dTraction = tNormalL * trans( tNormalL );
                    tContactTerm_dGap      = tNitscheParam * tNormalL;
                    tContactTerm_dNormal   = tNormalL * trans( tTraction ) - ( 2 * tPressure + tNitscheParam * tGap ) * tNormalL * trans( tNormalL ) + tContactTerm * tIdentity;
                }

                // compute the tangential projection operator (used to compute the derivative of the normal w.r.t the displacement)
                Matrix< DDRMat > const tTangentialProjector = tIdentity - tNormalL * trans( tNormalL );    // ( 2 x 2 )

                // compute the derivative of the normal w.r.t the variation of the displacement (Mlika (2018) Eq. (4.12))
                // TODO @ff: Can you interchange the variation of the displacement and the actual displacement?
                Matrix< DDRMat > const tNormal_dU    = -tTangentialProjector * trans( inv( tFL ) ) * trans( tGrad_UL ) * tNormalL;    // TODO @ff: Why is this matrix (2x1) and not (2x2)?
                Matrix< DDRMat > const tNormal_ddelU = -tTangentialProjector * trans( inv( tFL ) ) * tGrad_dUL * tNormalL;            // TODO @ff: tGrad_dUL is (2x4) but required is (2x2)

                // compute the derivative of the gap w.r.t the displacement
                const Matrix< DDRMat > tGap_dU = -( tNormalF / dot( tNormalF, tNormalL ) ) * ( tdUL - tdUF + tGap * tNormal_dU );    //

                // compute the derivative of the reference coordinate on the follower side (F) w.r.t the displacement variation (Mlika (2018) Eq. (4.13))
                const Matrix< DDRMat > tXF_dU = inv( tFF ) *                                                           //
                                                ( tIdentity - ( tNormalL * tNormalF ) / dot( tNormalL, tNormalF ) )    //
                                              * ( tdUL - tdUF - tGap * trans( inv( tFL ) ) * trans( tLeaderDofs->dnNdxn( 1 ) ) * tNormalL );

                // compute derivative with the contribution of the traction and the coordinate Y in the 4th line of Mlika (2018) Eq. (4.32)
                const Matrix< DDRMat > tTracYterm = mTheta * tTestTraction_dU + tNitscheParam * tLeaderDofs->dnNdxn( 1 ) * tXF_dU;

                auto tContribution =
                        -( mTheta / ( 2.0 * tNitscheParam ) ) * (                           //
                                trans( tTraction_dU ) * tNormalProjector * tTestTraction    // first integral in Mlika (2018) Eq. (4.32)
                                + tPressure * tTestTraction_dU                              // second integral in Mlika (2018) Eq. (4.32)
                                )
                        + ( 1.0 / ( 2.0 * tNitscheParam ) * (    //
                                    tContactTerm_dTraction * tTraction_dU + tContactTerm_dGap * tGap_dU + tContactTerm_dNormal * tNormal_dU )
                                        * ( mTheta * tTestTraction + tNitscheParam * ( tUF - tUL ) )    // third integral in Mlika (2018) Eq. (4.32)
                                + tContactTerm * tTracYterm                                             // fourth integral in Mlika (2018) Eq. (4.32)
                        );

                tJac( tResRange, tDepRange ) += aWStar * tContribution;


                // if ( tContactTerm < 0 )
                // {
                //     const Matrix< DDRMat > tdTractiondDof     = trans( mNormal ) * tConstitutiveModel->dTractiondDOF( tDofType, mNormal, CM_Function_Type::PK1 );
                //     const Matrix< DDRMat > tdTestTractiondDof = tConstitutiveModel->dTestTractiondDOF( tDofType, mNormal, tNormalProjector * tJump, tDisplDofTypes, CM_Function_Type::PK1 );
                //     tJac( tResRange, tDepRange ) += aWStar * ( -trans( tdUL ) * tNormalProjector * tdTractiondDof + mBeta * tdTestTractiondDof );
                // }
                // else
                // {
                //     const Matrix< DDRMat > tdTractiondDof     = tConstitutiveModel->dTractiondDOF( tDofType, mNormal, CM_Function_Type::PK1 );
                //     const Matrix< DDRMat > tdTestTractiondDof = tConstitutiveModel->dTestTractiondDOF( tDofType, mNormal, tNormalProjector * tTraction, tDisplDofTypes, CM_Function_Type::PK1 );
                //     tJac( tResRange, tDepRange ) += aWStar * -mBeta / tNitscheParam * ( trans( tTestTraction ) * tNormalProjector * tdTractiondDof + tdTestTractiondDof );
                // }
            }

            // if dependency of stabilization parameters on the dof type
            // if ( tStabilizationParameter->check_dof_dependency( tDofType, mtk::Leader_Follower::LEADER ) )
            // {
            //     // get the derivatives of the SPs
            //     const Matrix< DDRMat > tNitscheDer        = tStabilizationParameter->dSPdLeaderDOF( tDofType ).get_row( 0 );
            //     const Matrix< DDRMat > tLeaderWeightDer   = tStabilizationParameter->dSPdLeaderDOF( tDofType ).get_row( 1 );
            //     const Matrix< DDRMat > tFollowerWeightDer = tStabilizationParameter->dSPdLeaderDOF( tDofType ).get_row( 2 );
            //
            //     if ( tContactTerm < 0 )
            //     {
            //         // add contribution to Jacobian
            //         tJacMM += aWStar * (                                                                                                                            //
            //                           ( -trans(tdUL) * tNormalProjectorLeader * tTraction                                                                //
            //                                   + mBeta * trans(tTestTraction) * tNormalProjectorLeader * tJump    //
            //                                   + tNitsche * tLeaderDofs->N_trans() * tNormalProjectorLeader * tJump )                                          //
            //                                   * tLeaderWeightDer                                                                                                    //
            //                           + tLeaderDofs->N_trans() * tNormalProjectorLeader * tJump * tNitscheDer );
            //     }
            //     else
            //     {
            //         tJacMM += aWStar * (    //
            //                           -mBeta * tConstitutiveModel->testTraction_trans( mNormal, tDisplDofTypes ) * tNormalProjectorLeader * tTraction )
            //                 * ( 1.0 / tNitsche * tLeaderWeightDer - 1 / tNitsche / tNitsche * tNitscheDer );
            //     }
            // }
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian_and_residual( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian_and_residual - This function does nothing." );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_dRdp( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_dRdp - This function does nothing." );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
