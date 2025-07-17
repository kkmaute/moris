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
#include "fn_isfinite.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
#include "fn_trans.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include <memory>
#include <utility>

namespace moris::fem
{
    IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::IWG_Isotropic_Struc_Nonlinear_Contact_Mlika( sint aBeta )
            : mBeta( aBeta )    // sign for symmetric/unsymmetric Nitsche
            , mTheta( -mBeta )
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

        // test functions of leader (L) and follower (F)
        const Matrix< DDRMat > tNL = tLeaderDofs->N();      // ( 2 x 8 )
        const Matrix< DDRMat > tNF = tFollowerDofs->N();    // ( 2 x 8 )

        // coordinates of leader (L) and follower (F) in the reference configuration
        const Matrix< DDRMat > tXL = trans( tLeaderGeometry->valx() );      // ( 2 x 1 )
        const Matrix< DDRMat > tXF = trans( tFollowerGeometry->valx() );    // ( 2 x 1 )

        Matrix< DDRMat > const txL = tLeaderGeometry->valx_current( tLeaderDofs );
        Matrix< DDRMat > const txF = tFollowerGeometry->valx_current( tFollowerDofs );

        // compute the deformation gradients at the leader and follower side (all 2x1)
        const Matrix< DDRMat > tNormalCurrent   = tLeaderGeometry->get_normal_current( tLeaderDofs );
        const Matrix< DDRMat > tNormalReference = tLeaderGeometry->get_normal();
        const Matrix< DDRMat > tNormalProjector = tNormalCurrent * trans( tNormalCurrent );

        // gap between points in current configuration
        const Matrix< DDRMat > tJump = tUL - tUF - tXF + tXL;

        // evaluate traction and their pressure equivalent (normal component)
        const Matrix< DDRMat > tTraction = tConstitutiveModel->traction( tNormalReference, CM_Function_Type::PK1 );    // ( 2 x 1 )
        const real             tPressure = dot( tTraction, tNormalCurrent );

        // evaluate test traction and the corresponding pressure equivalent (normal component) by computing the dot product between
        // each column of tTestTraction and tNormalCurrent. This is done by element-wise multiplication and then summing the columns.
        const Matrix< DDRMat > tTestTraction = tConstitutiveModel->testTraction( tNormalReference, tDisplDofTypes, CM_Function_Type::PK1 );    // ( 2 x 8 )

        // Compute the contact term C_gamma(sigma, g, n), see Mlika (2018) Eq. (4.2.1). If the contact term is negative, the second term of the residual is not zero.
        const real tContactTerm = tPressure - tNitscheParam * dot( tJump, tNormalCurrent );
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

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                "IWG_Isotropic_Struc_Linear_Contact_Gap::compute_residual - Residual contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian( real aWStar )
    {
        MORIS_ERROR( false, "IWG_Isotropic_Struc_Nonlinear_Contact_Mlika::compute_jacobian - Jacobians for Nonconformal Contact can onl." );
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
