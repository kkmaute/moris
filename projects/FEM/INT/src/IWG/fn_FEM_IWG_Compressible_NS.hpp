/*
 * fn_FEM_IWG_Compressible_NS.hpp
 *
 *  Created on: Mar 18, 2021
 *      Author: wunsch
 * 
 * Free functions used in the IWGs for the unified compressible flow equations. 
 * Contains functions to compute various flux matrices, 
 * their derivatives, and functions to check Dof Lists
 */

#ifndef SRC_FEM_FN_FEM_IWG_COMPRESSIBLE_NS_HPP_
#define SRC_FEM_FN_FEM_IWG_COMPRESSIBLE_NS_HPP_

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        bool check_residual_dof_types( const moris::Cell< MSI::Dof_Type > & aResidualDofTypes );

        //------------------------------------------------------------------------------

        bool check_dof_dependencies(
                fem::Set                                           * aFemSet, 
                const moris::Cell< MSI::Dof_Type >                 & aResidualDofTypes,
                const moris::Cell< moris::Cell< MSI::Dof_Type > >  & aRequestedDofTypeList );   

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void eval_A( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                moris::Cell< Matrix< DDRMat > >       & aAMats );

        //------------------------------------------------------------------------------

        void eval_A0_DOF( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                moris::Cell< Matrix< DDRMat > >       & adA2dDOF );

        void eval_A1_DOF( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                moris::Cell< Matrix< DDRMat > >       & adA2dDOF );

        void eval_A2_DOF( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                moris::Cell< Matrix< DDRMat > >       & adA2dDOF );

        void eval_A3_DOF( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                moris::Cell< Matrix< DDRMat > >       & adA3dDOF );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void eval_K( 
                std::shared_ptr< Property >                      aPropDynamicViscosity,  
                std::shared_ptr< Property >                      aPropThermalConductivity,
                Field_Interpolator_Manager                     * aMasterFIManager,
                moris::Cell< moris::Cell< Matrix< DDRMat > > > & aKMats );

        void eval_KijYj( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                Matrix< DDRMat >                      & aKijYj );

        void eval_KijYji( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                Matrix< DDRMat >                      & aKijYji );

        //------------------------------------------------------------------------------

        void eval_KijYjDOF( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aDofType,
                moris::Cell< Matrix< DDRMat > >       & aKijYjDOF );

        void eval_KijYjiDOF( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aDofType,
                Matrix< DDRMat >                      & aKijYji );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        Matrix< DDRMat > unfold_flat_tensor( const Matrix< DDRMat > & aFlattenedTensor );

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_IWG_COMPRESSIBLE_NS_HPP_ */