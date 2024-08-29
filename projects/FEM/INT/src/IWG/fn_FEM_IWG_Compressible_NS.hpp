/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Compressible_NS.hpp
 *
 */

#ifndef SRC_FEM_FN_FEM_IWG_COMPRESSIBLE_NS_HPP_
#define SRC_FEM_FN_FEM_IWG_COMPRESSIBLE_NS_HPP_

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    bool check_residual_dof_types(
            const Vector< Vector< MSI::Dof_Type > > &aResidualDofTypes );

    //------------------------------------------------------------------------------

    bool check_dof_dependencies(
            fem::Set                                *aFemSet,
            const Vector< Vector< MSI::Dof_Type > > &aResidualDofTypes,
            const Vector< Vector< MSI::Dof_Type > > &aRequestedDofTypeList );

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void eval_A(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &aAMats );

    //------------------------------------------------------------------------------

    void eval_dAdY(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const uint                                   aAind,
            const uint                                   aYind,
            Matrix< DDRMat >                            &adAdY );

    void eval_VL_dAdY(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const Matrix< DDRMat >                      &aVL,
            const uint                                   aI,
            Matrix< DDRMat >                            &aVLdAdY );

    void eval_dAdY_VR(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const Matrix< DDRMat >                      &aVR,
            const uint                                   aI,
            Matrix< DDRMat >                            &adAdYVR );

    //------------------------------------------------------------------------------

    void eval_A0_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA0dDOF );

    void eval_A1_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA1dDOF );

    void eval_A2_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA2dDOF );

    void eval_A3_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA3dDOF );

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    void eval_K(
            const std::shared_ptr< Property >    &aPropDynamicViscosity,
            const std::shared_ptr< Property >    &aPropThermalConductivity,
            Field_Interpolator_Manager           *aLeaderFIManager,
            Vector< Vector< Matrix< DDRMat > > > &aK );

    void eval_dKijdxi(
            const std::shared_ptr< Property > &aPropDynamicViscosity,
            const std::shared_ptr< Property > &aPropThermalConductivity,
            Field_Interpolator_Manager        *aLeaderFIManager,
            Vector< Matrix< DDRMat > >        &adKijdxi );

    void eval_KijYj(
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Matrix< DDRMat >                            &aKijYj );

    void eval_KijYji(
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Matrix< DDRMat >                            &aKijYji );

    //------------------------------------------------------------------------------

    void eval_dKdY(
            const std::shared_ptr< Property >    &aPropDynamicViscosity,
            const std::shared_ptr< Property >    &aPropThermalConductivity,
            Field_Interpolator_Manager           *aLeaderFIManager,
            const uint                            aYind,
            Vector< Vector< Matrix< DDRMat > > > &adKdY );

    void eval_VL_dKdY(
            const std::shared_ptr< Property > &aPropDynamicViscosity,
            const std::shared_ptr< Property > &aPropThermalConductivity,
            Field_Interpolator_Manager        *aLeaderFIManager,
            const Matrix< DDRMat >            &aVL,
            const uint                         aI,
            const uint                         aJ,
            Matrix< DDRMat >                  &aVLdKdY );

    void eval_dKdY_VR(
            const std::shared_ptr< Property > &aPropDynamicViscosity,
            const std::shared_ptr< Property > &aPropThermalConductivity,
            Field_Interpolator_Manager        *aLeaderFIManager,
            const Matrix< DDRMat >            &aVR,
            const uint                         aI,
            const uint                         aJ,
            Matrix< DDRMat >                  &adKdYVR );

    void eval_VL_dKijidY(
            const std::shared_ptr< Property > &aPropDynamicViscosity,
            const std::shared_ptr< Property > &aPropThermalConductivity,
            Field_Interpolator_Manager        *aLeaderFIManager,
            const Matrix< DDRMat >            &aVR,
            const uint                         aJ,
            Vector< Matrix< DDRMat > >        &aVLdKijidY );

    void eval_dKijidY_VR(
            const std::shared_ptr< Property > &aPropDynamicViscosity,
            const std::shared_ptr< Property > &aPropThermalConductivity,
            Field_Interpolator_Manager        *aLeaderFIManager,
            const Matrix< DDRMat >            &aVR,
            const uint                         aJ,
            Vector< Matrix< DDRMat > >        &adKijidYVR );

    //------------------------------------------------------------------------------

    void eval_KijYjDOF(
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const Vector< MSI::Dof_Type >               &aDofType,
            Vector< Matrix< DDRMat > >                  &aKijYjDOF );

    void eval_KijYjiDOF(
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const Vector< MSI::Dof_Type >               &aDofType,
            Matrix< DDRMat >                            &aKijYjiDOF );

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    Matrix< DDRMat > unfold_flat_tensor( const Matrix< DDRMat > &aFlattenedTensor );

    //------------------------------------------------------------------------------

    uint convert_index_pair_to_flat( const uint aI, const uint aJ, const uint aNumSpaceDims );

    //------------------------------------------------------------------------------

}    // namespace moris::fem

#endif /* SRC_FEM_FN_FEM_IWG_COMPRESSIBLE_NS_HPP_ */
