/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IQI_HPP_
#define SRC_FEM_CL_FEM_IQI_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src
#include "cl_Matrix.hpp"         //LNA/src
// MRS/COR/src           // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
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
// FEM/VIS/src
#include "cl_VIS_Output_Enums.hpp"
// GEN/src
#include "GEN_Data_Types.hpp"
// LINALG/src
#include "fn_vectorize.hpp"
#include "cl_FEM_EvaluableTerm.hpp"

namespace moris
{
    namespace fem
    {
        class Set;
        class Cluster;
        class Field_Interpolator_Manager;

        /**
         * Integrand of a quantity of interest
         */
        class IQI : public EvaluableTerm
        {
          protected:
            enum fem::IQI_Type      mFEMIQIType;
            sint                    mIQITypeIndex = -1;
            const real              mToleranceFD  = 1e-12;    // tolerance for FD perturbation
            Vector< MSI::Dof_Type > mQuantityDofType;         // quantity dof type for IQI dof, max dof

            // local string to dof enum map
            std::map< std::string, MSI::Dof_Type > mLeaderDofMap;
            std::map< std::string, MSI::Dof_Type > mFollowerDofMap;

          private:
            // Normalization
            real        mReferenceValue    = 1.0;
            bool        mNormalized        = false;
            std::string mNormalizationType = "none";

            // function pointers
            void ( IQI::*m_compute_dQIdu_FD )(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType ) = nullptr;
            void ( IQI::*m_compute_dQIdp_FD_material )(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType ) = nullptr;
            void ( IQI::*m_compute_dQIdp_FD_geometry )(
                    real                          aWStar,
                    real                          aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices ) = nullptr;

            // function pointer for building the perturbation size for FD
            real ( IQI::*m_build_perturbation_size )(
                    const real& aPerturbation,
                    const real& aCoefficientToPerturb,
                    const real& aMaxPerturbation,
                    const real& aTolerance ) = nullptr;


          public:
            IQI()          = default;
            virtual ~IQI() = default;

            /**
             * get fem IQI type
             */
            enum fem::IQI_Type
            get_fem_IQI_type() { return mFEMIQIType; }

            /**
             * set quantity dof type (IQI dof, max dof)
             * @param[ in ] aQuantityDofType a cell of residual dof types
             */
            void
            set_quantity_dof_type( const Vector< MSI::Dof_Type >& aQuantityDofType ) { mQuantityDofType = aQuantityDofType; }

            /**
             * return a dof type for the quantity (IQI dof, max dof)
             * @param[ out ] mQuantityDofType a cell of residual dof types
             */
            const Vector< MSI::Dof_Type >&
            get_quantity_dof_type() const { return mQuantityDofType; }

            /**
             * Sets the reference values for this IQI.
             *
             * @param aReferenceValue Reference value for scaling the IQI, can be a norm if IQI is a vector.
             */
            void set_reference_value( real aReferenceValue );

            void set_normalization_type( std::string aNormalizationType ) { mNormalizationType = aNormalizationType; }

            std::string get_normalization_type() const { return mNormalizationType; }

            /*
             * set fem set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_fem_set( Set* aSet ) override
            {
                EvaluableTerm::set_fem_set( aSet );
                this->set_function_pointers();    // set function pointer for dQIdu, dQIdp
            }

            /*
             * set fem set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_function_pointers();

            /*
             * set output type index
             * @param[ in ] aOutputTypeIndex output type index
             */
            void
            set_output_type_index( sint aOutputTypeIndex )
            {
                mIQITypeIndex = aOutputTypeIndex;
            }

            /**
             * evaluate the derivative of the quantity of interest
             * wrt to requested dof types by finite difference
             * @param[ in ] aWStar        weight associated to the evaluation point
             * @param[ in ] aPerturbation real for relative perturbation of the dof values
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void
            compute_dQIdu_FD(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL )
            {
                // compute dQIdu geometry by FD
                ( this->*m_compute_dQIdu_FD )( aWStar, aPerturbation, aFDSchemeType );
            }

            void select_dQIdu_FD(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            /**
             * check the derivative of the quantity of interest wrt to dof types
             * with evaluation by finite difference
             * @param[ in ] aPerturbation real for perturbation of the dof values
             * @param[ in ] aEpsilon      real for tolerance
             * @param[ in ] adQIdu        matrix to fill with derivative of QI wrt dof types
             * @param[ in ] adQIduFD      matrix to fill with derivative of QI wrt dof types
             *                            evaluated by finite difference
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            bool check_dQIdu_FD(
                    real               aWStar,
                    real               aPerturbation,
                    real               aEpsilon,
                    Matrix< DDRMat >&  adQIdu,
                    Matrix< DDRMat >&  adQIduFD,
                    bool               aErrorPrint,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

            /**
             * evaluate the derivative of the quantity of interest wrt to dv types
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void
            compute_dQIdp( real aWStar )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdp - Not implemented for base class. " );
            }

            /**
             * evaluate the derivative of the quantity of interest
             * wrt to material dv by finite difference
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation dv relative perturbation
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void
            compute_dQIdp_FD_material(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL )
            {
                // compute dQIdp geometry by FD
                ( this->*m_compute_dQIdp_FD_material )( aWStar, aPerturbation, aFDSchemeType );
            }

            void select_dQIdp_FD_material(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            void
            select_dQIdp_FD_material_double(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType )
            {
                MORIS_ERROR( false, "IQI::select_dQIdp_FD_material_double - not implemented yet" );
            }

            /**
             * evaluate the derivative of the quantity of interest
             * wrt to geometry dv by finite difference
             * @param[ in ] aWStar            weight associated to evaluation point
             * @param[ in ] aPerturbation     pdv relative perturbation size
             * @param[ in ] aGeoLocalAssembly matrix filled with pdv local assembly indices
             * @param[ in ] aFDSchemeType     enum for FD scheme
             */
            void
            compute_dQIdp_FD_geometry(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices )
            {
                // compute dQIdp geometry by FD
                ( this->*m_compute_dQIdp_FD_geometry )(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType,
                        aGeoLocalAssembly,
                        aVertexIndices );
            }

            void select_dQIdp_FD_geometry_bulk(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices );

            void select_dQIdp_FD_geometry_sideset(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices );

            void
            select_dQIdp_FD_geometry_double(
                    moris::real                   aWStar,
                    moris::real                   aPerturbation,
                    fem::FDScheme_Type            aFDSchemeType,
                    Matrix< DDSMat >&             aGeoLocalAssembly,
                    Vector< Matrix< IndexMat > >& aVertexIndices )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdp_FD_geometry_double - not implemented yet" );
            }

            /**
             * add the contribution of the cluster measure derivatives to the derivative of
             * the quantity of interest wrt to geometry dv by finite difference
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation pdv relative perturbation size
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void add_cluster_measure_dQIdp_FD_geometry(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType );

            /**
             * Evaluate the quantity of interest.
             * @param[ out ] aQIVal quantity of interest matrix to fill
             */
            virtual void compute_QI( Matrix< DDRMat >& aQIVal ) = 0;

            /**
             * Evaluate the quantity of interest.
             * @param[ in ] aWStar            weight associated to evaluation point
             */
            virtual void compute_QI( real aWStar ) = 0;

            /**
             * Compute the derivative of the quantities of interest wrt requested dof types.
             * @param[ in ]  aDofType Dof type being evaluated
             * @param[ out ] adQIdu derivative of quantity of interest
             */
            virtual void compute_dQIdu(
                    Vector< MSI::Dof_Type >& aDofType,
                    Matrix< DDRMat >&        adQIdu ) = 0;

            /**
             * Compute the derivative of the quantities of interest
             * @param[ in ] aWStar weight associated to evaluation point
             */
            virtual void compute_dQIdu( real aWstar ) = 0;

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
             * @param[ out ] tDeltaH              perturbation size built for finite difference
             */
            real check_ig_coordinates_inside_ip_element(
                    const real&         aPerturbation,
                    const real&         aCoefficientToPerturb,
                    const uint&         aSpatialDirection,
                    fem::FDScheme_Type& aUsedFDScheme );

            /**
             * get the matrix dimension of the IQI in order to initialize the size of the mGloblaIQIVal
             * returns 1*1 which is a scaler by default
             */
            virtual std::pair< uint, uint >
            get_matrix_dim()
            {
                return std::make_pair( 1, 1 );
            }
        };
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IQI_HPP_ */
