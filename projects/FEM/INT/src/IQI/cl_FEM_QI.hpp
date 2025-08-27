/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_QI.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_QI_HPP_
#define SRC_FEM_CL_FEM_QI_HPP_

#include <utility>

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

namespace moris::fem
{
    class Set;
    class Cluster;
    class Field_Interpolator_Manager;

    //------------------------------------------------------------------------------
    /**
     * Integrand of a quantity of interest
     */
    class QI
    {
      protected:
        // FEM QI type
        enum fem::IQI_Type mFEMQIType;

        // QI type index
        sint mQITypeIndex = -1;

        // normal
        Matrix< DDRMat > mNormal;

        // tolerance for FD perturbation
        const real mToleranceFD = 1e-12;

        // flag for building global dv type list
        bool mGlobalDvBuild = true;

        bool mGlobalFieldBuild = true;

        // parameters
        Vector< Matrix< DDRMat > > mParameters;

        // QI name
        std::string mName;

      private:
        // Normalization
        real mReferenceValue = 1.0;
        bool mNormalized     = false;

        // function pointers
        void ( QI::*m_compute_dQIdu_FD )(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType ) = nullptr;
        void ( QI::*m_compute_dQIdp_FD_geometry )(
                real                          aWStar,
                real                          aPerturbation,
                fem::FDScheme_Type            aFDSchemeType,
                Matrix< DDSMat >&             aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices ) = nullptr;

        // function pointer for building the perturbation size for FD
        real ( QI::*m_build_perturbation_size )(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance ) = nullptr;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * constructor
         */
        QI() {};

        //------------------------------------------------------------------------------
        /**
         * destructor
         */
        virtual ~QI() {};

        //------------------------------------------------------------------------------
        /**
         * set name
         * param[ in ] aName a string for CM name
         */
        void
        set_name( std::string aName )
        {
            mName = std::move( aName );
        }

        //------------------------------------------------------------------------------
        /**
         * get name
         * param[ out ] mName a string for CM name
         */
        std::string
        get_name()
        {
            return mName;
        }

        //------------------------------------------------------------------------------
        /**
         * print names
         */
        void print_names();

        //------------------------------------------------------------------------------
        /**
         * get fem QI type
         */
        enum fem::IQI_Type
        get_fem_QI_type()
        {
            return mFEMQIType;
        }

        //------------------------------------------------------------------------------
        /**
         * rest evaluation flags for the QI
         */
        void reset_eval_flags();

        //------------------------------------------------------------------------------

        /**
         * Sets the reference values for this QI.
         *
         * @param aReferenceValue Reference value for scaling the QI, can be a norm if QI is a vector.
         */
        void set_reference_value( real aReferenceValue );

        //------------------------------------------------------------------------------
        /*
         * set fem set pointer
         * @param[ in ] aSetPointer a FEM set pointer
         */
        void set_function_pointers();

        //------------------------------------------------------------------------------
        /*
         * set output type index
         * @param[ in ] aOutputTypeIndex output type index
         */
        void
        set_output_type_index( sint aOutputTypeIndex )
        {
            mQITypeIndex = aOutputTypeIndex;
        }

        //------------------------------------------------------------------------------
        /**
         * set normal
         * @param[ in ] aNormal normal vector
         */
        void
        set_normal( const Matrix< DDRMat >& aNormal )
        {
            mNormal = aNormal;
        }

        //------------------------------------------------------------------------------
        /**
         * set parameters
         * @param[ in ] aParameters a list of parameters
         */
        virtual void
        set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
        {
            // set a cluster
            mParameters = aParameters;
        }

        //------------------------------------------------------------------------------
        /*
         * set field interpolator manager pointer
         * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
         */
        void set_field_interpolator_manager(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * set field interpolator manager for previous time step
         * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
         * @param[ in ] aIsLeader                 an enum for leader or follower
         */
        void set_field_interpolator_manager_previous_time(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * set field interpolator manager for eigen vector
         * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
         * @param[ in ] aIsLeader                 an enum for leader or follower
         */
        void set_field_interpolator_manager_eigen_vector(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /*
         * get field interpolator manager
         * @param[ out ] aFieldInterpolatorManager a field interpolator manager pointer
         * @param[ in ]  aIsLeader                 an enum for leader or follower
         */
        Field_Interpolator_Manager* get_field_interpolator_manager(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set dof type list for QI
         * @param[ in ] aDofTypes a cell of group of dof types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                mtk::Leader_Follower                     aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set dof types
         * @param[ in ] aDofTypes   list of group of dof types
         * @param[ in ] aDofStrings list of strings describing the group dof types
         * @param[ in ] aIsLeader   enum for leader or follower
         */
        virtual void
        set_dof_type_list(
                Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                Vector< std::string >&             aDofStrings,
                mtk::Leader_Follower               aIsLeader = mtk::Leader_Follower::LEADER )
        {
            MORIS_ERROR( false, "QI::set_dof_type_list - not implemented for base class." );
        }

        //------------------------------------------------------------------------------
        /**
         * return a cell of dof types
         * @param[ in ] aIsLeader enum leader or follower
         */
        const Vector< Vector< MSI::Dof_Type > >& get_dof_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get a non unique list of dof type including
         * QI, property, constitutive and stabilization dependencies
         * for both leader and follower
         */
        void get_non_unique_dof_dv_and_field_types(
                Vector< Vector< MSI::Dof_Type > >&   aDofTypes,
                Vector< Vector< gen::PDV_Type > >&   aDvTypes,
                Vector< Vector< mtk::Field_Type > >& aFieldTypes );

        //------------------------------------------------------------------------------
        /**
         * get a unique global dof type list including
         * QI, property, constitutive and stabilization dependencies
         * @param[ in ] aIsLeader enum leader or follower
         */
        const Vector< Vector< MSI::Dof_Type > >& get_global_dof_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * build a global dof and dv type lists including
         * QI, property, constitutive and stabilization dependencies
         * ( a list for leader and a list for follower dof and dv types )
         */
        void build_global_dof_dv_and_field_type_list();

        //------------------------------------------------------------------------------
        /**
         * set dv types
         * @param[ in ] aDvTypes  a cell of cell of dv types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_dv_type_list(
                const Vector< Vector< gen::PDV_Type > >& aDvTypes,
                mtk::Leader_Follower                     aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * return a cell of dv types
         * @param[ in ] aIsLeader enum leader or follower
         */
        Vector< Vector< gen::PDV_Type > >& get_dv_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * return a cell of field types active for the QI
         * @param[ in ]  aIsLeader enum leader or follower
         * @param[ out ] aFieldTypes a list of group of field types
         */
        const Vector< Vector< mtk::Field_Type > >& get_field_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //------------------------------------------------------------------------------
        /**
         * set QI active field types
         * @param[ in ] aFieldTypes a list of group of field types
         * @param[ in ] aIsLeader enum for leader or follower
         */
        void set_field_type_list(
                const Vector< Vector< mtk::Field_Type > >& aDvTypes,
                mtk::Leader_Follower                       aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get a non unique list of dv type including
         * QI, property, constitutive and stabilization dependencies
         * @param[ in ] aGlobalDvTypeList a non unique list of dv types to fill
         */
        void get_non_unique_global_dv_type_list( Vector< gen::PDV_Type >& aGlobalDvTypeList );

        //------------------------------------------------------------------------------
        /**
         * get a unique global dv type list including
         * QI, property, constitutive and stabilization dependencies
         * @param[ in ] aIsLeader enum leader or follower
         */
        const Vector< Vector< gen::PDV_Type > >& get_global_dv_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get global field type list. TODO: Field types are only used by the IWG.
         * If a user wants to use them in a property or CM this cuntion has to be modified in the same way than get_global_dof_type_list()
         * @param[ in ]  aIsLeader    enum leader or follower
         * @param[ out ] mFieldTypes global list of group of dv types
         */
        const Vector< Vector< mtk::Field_Type > >& get_global_field_type_list(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set property
         * @param[ in ] aProperty       a property pointer
         * @param[ in ] aPropertyString a string describing the property
         * @param[ in ] aIsLeader       enum leader or follower
         */
        void set_property(
                std::shared_ptr< Property > aProperty,
                const std::string&          aPropertyString,
                mtk::Leader_Follower        aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get leader or follower properties
         * @param[ in ]  aIsLeader   enum leader or follower
         * @param[ out ] aProperties cell of property pointers
         */
        Vector< std::shared_ptr< fem::Property > >& get_properties(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set constitutive model
         * @param[ in ] aConstitutiveModel  a constitutive model pointer
         * @param[ in ] aConstitutiveString a string describing the constitutive model
         * @param[ in ] aIsLeader           an enum for leader or follower
         */
        void set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                const std::string&                    aConstitutiveString,
                mtk::Leader_Follower                  aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * get leader or follower constitutive models
         * @param[ in ]  aIsLeader           enum leader or follower
         * @param[ out ] aConstitutiveModels cell of constitutive model pointers
         */
        Vector< std::shared_ptr< fem::Constitutive_Model > >& get_constitutive_models(
                mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------
        /**
         * set stabilization parameter
         * @param[ in ] aStabilizationParameter a stabilization parameter pointer
         * @param[ in ] aStabilizationString    a string defining the stabilization parameter
         */
        void set_stabilization_parameter(
                const std::shared_ptr< Stabilization_Parameter >& aStabilizationParameter,
                const std::string&                                aStabilizationString );

        //------------------------------------------------------------------------------
        /**
         * build requested dof type list
         */
        void build_requested_dof_type_list();

        //------------------------------------------------------------------------------
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

        //------------------------------------------------------------------------------
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

        //------------------------------------------------------------------------------
        /**
         * evaluate the derivative of the quantity of interest wrt to dv types
         * @param[ in ] aWStar weight associated to the evaluation point
         */
        virtual void
        compute_dQIdp( real aWStar )
        {
            MORIS_ERROR( false, "QI::compute_dQIdp - Not implemented for base class. " );
        }

        //------------------------------------------------------------------------------

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
            MORIS_ERROR( false, "QI::select_dQIdp_FD_material_double - not implemented yet" );
        }

        //------------------------------------------------------------------------------
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
            MORIS_ERROR( false, "QI::compute_dQIdp_FD_geometry_double - not implemented yet" );
        }

        //--------------------------------------------------f----------------------------
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

        //------------------------------------------------------------------------------
        /**
         * Evaluate the quantity of interest.
         * @param[ out ] aQIVal quantity of interest matrix to fill
         */
        virtual void compute_QI( Matrix< DDRMat >& aQIVal ) = 0;

        //------------------------------------------------------------------------------
        /**
         * Compute the derivative of the quantities of interest wrt requested dof types.
         * @param[ in ]  aDofType Dof type being evaluated
         * @param[ out ] adQIdu derivative of quantity of interest
         */
        virtual void compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu ) = 0;

        //------------------------------------------------------------------------------
        /**
         * Compute the derivative of the quantities of interest
         * @param[ in ] aWStar weight associated to evaluation point
         */
        virtual void compute_dQIdu( real aWstar ) = 0;

        //------------------------------------------------------------------------------
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

        //------------------------------------------------------------------------------
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

        //------------------------------------------------------------------------------
        /**
         * get the matrix dimension of the QI in order to initialize the size of the mGloblaQIVal
         * returns 1*1 which is a scaler by default
         */
        virtual std::pair< uint, uint >
        get_matrix_dim()
        {
            return std::make_pair( 1, 1 );
        }
    };
    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_QI_HPP_ */
