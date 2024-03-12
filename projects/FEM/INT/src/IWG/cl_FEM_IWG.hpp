/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_
// MRS/CNT/src
#include "cl_Vector.hpp"
// LNA/src
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "fn_vectorize.hpp"
#include "fn_isfinite.hpp"
// MRS/COR/src // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Enums.hpp"
#include "fn_FEM_FD_Scheme.hpp"
// FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
// GEN/src
#include "GEN_Data_Types.hpp"
#include "cl_FEM_EvaluableTerm.hpp"

namespace moris
{
    namespace fem
    {
        class Set;
        class Cluster;
        class Field_Interpolator_Manager;
        class FEM_Model;

        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG : public EvaluableTerm
        {
          protected:
            // nodal weak BCs
            Matrix< DDRMat > mNodalWeakBCs;


            // residual dof type
            Vector< Vector< MSI::Dof_Type > > mResidualDofType;

            // Local string to int map for material models

            // interpolation order for IWG
            uint mOrder = MORIS_UINT_MAX;

            // tolerance for FD perturbation
            const real mToleranceFD = 1e-12;

            // bool for ghost
            bool mIsGhost = false;

          protected:
            // compute the jacobian using finite differencing (independent on the setting for the whole FEM set)
            bool mIsFDJacobian = false;

          protected:
            //! string for IWG name
            enum moris::fem::IWG_Type mIWGType = moris::fem::IWG_Type::UNDEFINED;

            // function pointers
            void ( IWG::*m_compute_jacobian_FD )(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType,
                    bool               aUseAbsolutePerturbations ) = nullptr;
            void ( IWG::*m_compute_dRdp_FD_material )(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType ) = nullptr;
            void ( IWG::*m_compute_dRdp_FD_geometry )(
                    real                          aWStar,
                    real                          aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices ) = nullptr;

            // function pointer for building the perturbation size for FD
            real ( IWG::*m_build_perturbation_size )(
                    const real& aPerturbation,
                    const real& aCoefficientToPerturb,
                    const real& aMaxPerturbation,
                    const real& aTolerance ) = nullptr;


          public:
            IWG()          = default;
            virtual ~IWG() = default;

            /**
             * set ghost flag
             * param[ in ] aIsGhost bool true if IWG for ghost
             */
            void set_ghost_flag( bool aIsGhost ) { mIsGhost = aIsGhost; }

            /**
             * get ghost flag
             * param[ out ] mIsGhost bool true if IWG for ghost
             */
            bool
            get_ghost_flag() { return mIsGhost; }

            /**
             * get IWG type
             * param[ out ] mIWGType an enum of the IWG type. type only implemented for TIME_CONTINUITY_DOF. All others return UNDEFINED.
             *              If needed you can implement the type for the others. Just folow the TIME_CONTINUITY_DOF
             */
            enum moris::fem::IWG_Type
            get_IWG_type() { return mIWGType; };


            bool is_analytical_jacobian() const { return mIsAnalyticalJacobian; };
            void set_is_analytical_jacobian( bool aIsAnalyticalJacobian ) { mIsAnalyticalJacobian = aIsAnalyticalJacobian; };

            bool is_fd_jacobian() const
            {
                return mIsFDJacobian;
            };

            void set_is_fd_jacobian( bool aIsFDJacobian )
            {
                mIsFDJacobian = aIsFDJacobian;
            };
            //------------------------------------------------------------------------------
            /*
             * set member set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_fem_set( Set* aSet ) override
            {
                EvaluableTerm::set_fem_set( aSet );
                this->set_function_pointers();    // set function pointer for dQIdu, dQIdp
            }

            /*
             * set fem cluster pointer
             * @param[ in ] aClusterPointer a FEM cluster pointer
             */
            void
            set_cluster_pointer( fem::Cluster* aClusterPointer ) { mCluster = aClusterPointer; }

            /*
             * set function pointers
             */
            void set_function_pointers();

            /*
             * get member set pointer
             * @param[ out ] aSetPointer a FEM set pointer
             */
            Set*
            get_set_pointer() { return mSet; }

            /*
             * free memory
             */
            void
            free_memory() {}

            /**
             * set nodal weak BCs
             * @param[ in ] aNodalWeakBCs matrix with nodal values
             */
            void
            set_nodal_weak_bcs( const Matrix< DDRMat >& aNodalWeakBCs ) { mNodalWeakBCs = aNodalWeakBCs; }

            /**
             * set residual dof type
             * @param[ in ] aResidualdofType a cell of residual dof types
             */
            void
            set_residual_dof_type( const Vector< Vector< MSI::Dof_Type > >& aResidualDofType ) { mResidualDofType = aResidualDofType; }

            /**
             * return a dof type for the residual
             * @param[ out ] aResidualdofType a cell of residual dof types
             */
            const Vector< Vector< MSI::Dof_Type > >&
            get_residual_dof_type() const { return mResidualDofType; }

            /**
             * set interpolation order for the residual dof type
             */
            void set_interpolation_order();
            void set_interpolation_order( uint aOrder ) { mOrder = aOrder; }

            /**
             * check that field interpolators were assigned
             * @param[ in ]  aIsLeader enum leader or follower
             */
            void check_field_interpolators( mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

            /**
             * evaluate the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_residual( real aWStar ) = 0;

            /**
             * evaluate the Jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_jacobian( real aWStar ) = 0;

            virtual void
            compute_jacobian_previous( real aWStar )
            {
                MORIS_ERROR( false, "compute_jacobian_previous() not implemented" );
            };

            /**
             * evaluate the Jacobian by finite difference
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void
            compute_jacobian_FD(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType             = fem::FDScheme_Type::POINT_5,
                    bool               aUseAbsolutePerturbations = false )
            {
                // compute jacobian by FD
                ( this->*m_compute_jacobian_FD )( aWStar, aPerturbation, aFDSchemeType, aUseAbsolutePerturbations );
            }

            void select_jacobian_FD(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType,
                    bool               aUseAbsolutePerturbations );

            void select_jacobian_FD_double(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType,
                    bool               aUseAbsolutePerturbations );

            /**
             * evaluate the residual and the Jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_jacobian_and_residual( real aWStar ) = 0;

            /**
             * check the Jacobian with FD
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aEpsilon      real for check
             * @param[ in ] aWStar        real weight associated to evaluation point
             * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
             * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
             * @param[ in ] aErrorPrint   bool set to true to print non matching values in jacobian
             */
            bool check_jacobian(
                    real              aPerturbation,
                    real              aEpsilon,
                    real              aWStar,
                    Matrix< DDRMat >& aJacobians,
                    Matrix< DDRMat >& aJacobiansFD,
                    bool              aErrorPrint               = false,
                    bool              aUseAbsolutePerturbations = false );

            /**
             * check Jacobian that uses multiple dof types with FD
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aEpsilon      real for check
             * @param[ in ] aWStar        real weight associated to evaluation point
             * @param[ in ] aJacobians    cell of cell of matrices to fill with Jacobians
             * @param[ in ] aJacobians_FD cell of cell of matrices to fill with Jacobians by FD
             * @param[ in ] aErrorPrint   bool set to true to print non matching values in jacobian
             */
            bool check_jacobian_multi_residual(
                    real              aPerturbation,
                    real              aEpsilon,
                    real              aWStar,
                    Matrix< DDRMat >& aJacobians,
                    Matrix< DDRMat >& aJacobiansFD,
                    bool              aErrorPrint               = false,
                    bool              aMaxErrorPrint            = false,
                    moris::real       aFDtolerance              = -1.0,
                    bool              aUseAbsolutePerturbations = false );

            /**
             * evaluate the derivative of the residual wrt the design variables
             * @param[ in ] aWStar weight associated to evaluation point
             */
            virtual void compute_dRdp( real aWStar ) = 0;

            /**
             * evaluate the derivative of the residual
             * wrt the material design variables by finite difference
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation real for dv perturbation
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void
            compute_dRdp_FD_material(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL )
            {
                // compute jacobian by FD
                ( this->*m_compute_dRdp_FD_material )( aWStar, aPerturbation, aFDSchemeType );
            }

            void select_dRdp_FD_material(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            void select_dRdp_FD_material_double(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            /**
             * evaluate the derivative of the residual
             * wrt the geometry design variables by finite difference
             * @param[ in ] aWStar            weight associated to evaluation point
             * @param[ in ] aPerturbation     real for relative dv perturbation
             * @param[ in ] aGeoLocalAssembly matrix filled with pdv local assembly indices
             * @param[ in ] aFDSchemeType     enum for FD scheme
             */
            void
            compute_dRdp_FD_geometry(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices )
            {
                // compute jacobian by FD
                ( this->*m_compute_dRdp_FD_geometry )(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType,
                        aGeoLocalAssembly,
                        aVertexIndices );
            }

            void select_dRdp_FD_geometry_bulk(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices );

            void select_dRdp_FD_geometry_sideset(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices );

            void select_dRdp_FD_geometry_time_sideset(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices );

            void select_dRdp_FD_geometry_double(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices );

            /**
             * add the contribution of the cluster measure derivatives to the derivative of
             * the quantity of interest wrt to geometry dv by finite difference
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation pdv relative perturbation size
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void add_cluster_measure_dRdp_FD_geometry(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            /**
             * add the contribution of the cluster measure derivatives to the derivative of
             * the quantity of interest wrt to geometry dv by finite difference
             * for double sideset
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation pdv relative perturbation size
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void add_cluster_measure_dRdp_FD_geometry_double(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            /**
             * build perturbation size for finite difference
             * @param[ in ] aPerturbation         provided perturbation size from input
             * @param[ in ] aCoefficientToPerturb coefficient to perturb
             * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
             */
            real build_perturbation_size(
                    const real& aPerturbation,
                    const real& aCoefficientToPerturb,
                    const real& aMaxPerturbation,
                    const real& aTolerance );

            /**
             * build relative perturbation size for finite difference
             * @param[ in ] aPerturbation         provided perturbation size from input
             * @param[ in ] aCoefficientToPerturb coefficient to perturb
             * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
             */
            real build_perturbation_size_relative(
                    const real& aPerturbation,
                    const real& aCoefficientToPerturb,
                    const real& aMaxPerturbation,
                    const real& aTolerance );

            /**
             * build absolute perturbation size for finite difference
             * @param[ in ] aPerturbation         provided perturbation size from input
             * @param[ in ] aCoefficientToPerturb coefficient to perturb
             * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
             */
            real build_perturbation_size_absolute(
                    const real& aPerturbation,
                    const real& aCoefficientToPerturb,
                    const real& aMaxPerturbation,
                    const real& aTolerance );

            /**
             * check if ig node still inside ip element after perturbation in a specific
             * space direction, if not adapt the finite difference scheme used
             * @param[ in ] aPerturbation         provided perturbation size from input
             * @param[ in ] aCoefficientToPerturb coefficient to perturb
             * @param[ in ] aSpatialDirection     spatial direction in which we perturb
             * @param[ in ] aUsedFDScheme         FD scheme to be used, updated
             * @param[ out ] aDeltaH              perturbation size built for finite difference
             */
            real check_ig_coordinates_inside_ip_element(
                    const real&         aPerturbation,
                    const real&         aCoefficientToPerturb,
                    const uint&         aSpatialDirection,
                    fem::FDScheme_Type& aUsedFDScheme );

            /**
             * reset evaluation flags
             */
            void reset_eval_flags() override;

            /**
             * reset evaluation flags specific to child IWG
             */
            virtual void reset_spec_eval_flags() {};

            Matrix< DDRMat > get_deformed_node_coordinates( Geometry_Interpolator* aGeometryInterpolator, Field_Interpolator* aFieldInterpolator ) const;
            Matrix< DDRMat > remap_nonconformal_rays( Field_Interpolator* aLeaderFieldInterpolator, Field_Interpolator* aFollowerFieldInterpolator, bool aDebugOutput ) const;
        };

    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
