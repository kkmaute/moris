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

#include "typedefs.hpp"     //MRS/COR/src
#include "cl_Cell.hpp"      //MRS/CON/src
#include "cl_Matrix.hpp"    //LNA/src
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
#include "cl_GEN_Pdv_Enums.hpp"
// LINALG/src
#include "fn_vectorize.hpp"

namespace moris
{
    namespace fem
    {
        class Set;
        class Cluster;
        class Field_Interpolator_Manager;

        //------------------------------------------------------------------------------
        /**
         * Integrand of a quantity of interest
         */
        class IQI
        {
          protected:
            // FEM set pointer
            fem::Set* mSet = nullptr;

            // cluster pointer
            fem::Cluster* mCluster = nullptr;

            // FEM IQI type
            enum fem::IQI_Type mFEMIQIType;

            // IQI type index
            sint mIQITypeIndex = -1;

            // normal
            Matrix< DDRMat > mNormal;

            // tolerance for FD perturbation
            const real mToleranceFD = 1e-12;

            // quantity dof type for IQI dof, max dof
            moris::Cell< MSI::Dof_Type > mQuantityDofType;

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave requested global dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedSlaveGlobalDofTypes;

            // flag for building global dof type list
            bool mGlobalDofBuild = true;

            // field interpolator manager pointer
            Field_Interpolator_Manager* mMasterFIManager         = nullptr;
            Field_Interpolator_Manager* mSlaveFIManager          = nullptr;
            Field_Interpolator_Manager* mMasterPreviousFIManager = nullptr;

            // master and slave dv type lists
            moris::Cell< moris::Cell< PDV_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< PDV_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< PDV_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< PDV_Type > > mSlaveGlobalDvTypes;

            // flag for building global dv type list
            bool mGlobalDvBuild = true;

            // master and slave field type lists
            moris::Cell< moris::Cell< mtk::Field_Type > > mMasterFieldTypes;
            moris::Cell< moris::Cell< mtk::Field_Type > > mSlaveFieldTypes;

            moris::Cell< moris::Cell< mtk::Field_Type > > mMasterGlobalFieldTypes;
            moris::Cell< moris::Cell< mtk::Field_Type > > mSlaveGlobalFieldTypes;

            bool mGlobalFieldBuild = true;

            // master and slave properties
            moris::Cell< std::shared_ptr< fem::Property > > mMasterProp;
            moris::Cell< std::shared_ptr< fem::Property > > mSlaveProp;

            // local string to int map for properties
            std::map< std::string, uint > mPropertyMap;

            // master and slave constitutive models
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mMasterCM;
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mSlaveCM;

            // Local string to int map for constitutive models
            std::map< std::string, uint > mConstitutiveMap;

            // stabilization parameters
            moris::Cell< std::shared_ptr< Stabilization_Parameter > > mStabilizationParam;

            // local string to int map for stabilizations
            std::map< std::string, uint > mStabilizationMap;

            // active cluster measure on IQI flag
            bool mActiveCMEAFlag = false;

            // local string to dof enum map
            std::map< std::string, MSI::Dof_Type > mMasterDofMap;
            std::map< std::string, MSI::Dof_Type > mSlaveDofMap;

            // parameters
            moris::Cell< Matrix< DDRMat > > mParameters;

            // IQI name
            std::string mName;

            // bulk type
            fem::Element_Type mBulkType = fem::Element_Type::BULK;

            // strings for master and slave phase name
            std::string mMasterPhaseName;
            std::string mSlavePhaseName;

            // bool for time continuity
            bool mTimeContinuity = false;

            // bool for time boundary
            bool mTimeBoundary = false;

          private:
            // Normalization
            real mReferenceValue = 1.0;
            bool mNormalized     = false;

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
                    real                               aWStar,
                    real                               aPerturbation,
                    fem::FDScheme_Type                 aFDSchemeType,
                    Matrix< DDSMat >&                  aGeoLocalAssembly,
                    moris::Cell< Matrix< IndexMat > >& aVertexIndices ) = nullptr;

            // function pointer for building the perturbation size for FD
            real ( IQI::*m_build_perturbation_size )(
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
            IQI(){};

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~IQI(){};

            //------------------------------------------------------------------------------
            /**
             * set name
             * param[ in ] aName a string for CM name
             */
            void
            set_name( std::string aName )
            {
                mName = aName;
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
             * set time continuity flag
             * param[ in ] aTimeContinuity bool true if IWG for time continuity
             */
            void
            set_time_continuity( bool aTimeContinuity )
            {
                mTimeContinuity = aTimeContinuity;
            }

            //------------------------------------------------------------------------------
            /**
             * get time continuity flag
             * param[ out ] mTimeContinuity ool true if IWG for time continuity
             */
            bool
            get_time_continuity()
            {
                return mTimeContinuity;
            }

            //------------------------------------------------------------------------------
            /**
             * set time boundary flag
             * param[ in ] aTimeBoundary bool true if IWG for time boundary
             */
            void
            set_time_boundary( bool aTimeBoundary )
            {
                mTimeBoundary = aTimeBoundary;
            }

            //------------------------------------------------------------------------------
            /**
             * get time boundary flag
             * param[ out ] mTimeBoundary bool true if IWG for time boundary
             */
            bool
            get_time_boundary()
            {
                return mTimeBoundary;
            }

            //------------------------------------------------------------------------------
            /**
             * set bulk type
             * @param[ in ] aBulkType element type for the IWG
             */
            void
            set_bulk_type( fem::Element_Type aBulkType )
            {
                mBulkType = aBulkType;
            }

            //------------------------------------------------------------------------------
            /**
             * get bulk type
             * @param[ out ] mBulkType element type for the IWG
             */
            fem::Element_Type
            get_bulk_type()
            {
                return mBulkType;
            }

            //------------------------------------------------------------------------------
            /**
             * get fem IQI type
             */
            enum fem::IQI_Type
            get_fem_IQI_type()
            {
                return mFEMIQIType;
            }

            //------------------------------------------------------------------------------
            /**
             * set phase name
             * param[ in ] aPhaseName a string for phase name
             * param[ in ] aIsMaster  an enum for master or slave
             */
            void set_phase_name(
                    std::string       aPhaseName,
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get phase name
             * param[ in ]  aIsMaster an enum for master or slave
             * param[ out ] mName     a string for phase name
             */
            std::string get_phase_name(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set quantity dof type (IQI dof, max dof)
             * @param[ in ] aQuantityDofType a cell of residual dof types
             */
            void
            set_quantity_dof_type( const moris::Cell< MSI::Dof_Type >& aQuantityDofType )
            {
                mQuantityDofType = aQuantityDofType;
            }

            //------------------------------------------------------------------------------
            /**
             * return a dof type for the quantity (IQI dof, max dof)
             * @param[ out ] mQuantityDofType a cell of residual dof types
             */
            const moris::Cell< MSI::Dof_Type >&
            get_quantity_dof_type() const
            {
                return mQuantityDofType;
            }

            //------------------------------------------------------------------------------
            /**
             * rest evaluation flags for the IQI
             */
            void reset_eval_flags();

            //------------------------------------------------------------------------------

            /**
             * Sets the reference values for this IQI.
             *
             * @param aReferenceValue Reference value for scaling the IQI, can be a norm if IQI is a vector.
             */
            void set_reference_value( real aReferenceValue );

            //------------------------------------------------------------------------------
            /*
             * set fem set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void
            set_set_pointer( Set* aSetPointer )
            {
                mSet = aSetPointer;

                // set function pointer for dQIdu, dQIdp
                this->set_function_pointers();
            }

            //------------------------------------------------------------------------------
            /*
             * set fem cluster pointer
             * @param[ in ] aClusterPointer a FEM cluster pointer
             */
            void
            set_cluster_pointer( fem::Cluster* aClusterPointer )
            {
                mCluster = aClusterPointer;
            }

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
                mIQITypeIndex = aOutputTypeIndex;
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
            void
            set_parameters( const moris::Cell< Matrix< DDRMat > >& aParameters )
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
                    mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /*
             * set field interpolator manager for previous time step
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             * @param[ in ] aIsMaster                 an enum for master or slave
             */
            void set_field_interpolator_manager_previous_time(
                    Field_Interpolator_Manager* aFieldInterpolatorManager,
                    mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /*
             * get field interpolator manager
             * @param[ out ] aFieldInterpolatorManager a field interpolator manager pointer
             * @param[ in ]  aIsMaster                 an enum for master or slave
             */
            Field_Interpolator_Manager* get_field_interpolator_manager(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set dof type list for IQI
             * @param[ in ] aDofTypes a cell of group of dof types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dof_type_list(
                    const moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                    mtk::Master_Slave                                  aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes   list of group of dof types
             * @param[ in ] aDofStrings list of strings describing the group dof types
             * @param[ in ] aIsMaster   enum for master or slave
             */
            virtual void
            set_dof_type_list(
                    moris::Cell< moris::Cell< MSI::Dof_Type > >& aDofTypes,
                    moris::Cell< std::string >&                  aDofStrings,
                    mtk::Master_Slave                            aIsMaster = mtk::Master_Slave::MASTER )
            {
                MORIS_ERROR( false, "IQI::set_dof_type_list - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] aIsMaster enum master or slave
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > >& get_dof_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get a non unique list of dof type including
             * IQI, property, constitutive and stabilization dependencies
             * for both master and slave
             */
            void get_non_unique_dof_dv_and_field_types(
                    moris::Cell< moris::Cell< MSI::Dof_Type > >&   aDofTypes,
                    moris::Cell< moris::Cell< PDV_Type > >&        aDvTypes,
                    moris::Cell< moris::Cell< mtk::Field_Type > >& aFieldTypes );

            //------------------------------------------------------------------------------
            /**
             * get a unique global dof type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > >& get_global_dof_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * build a global dof and dv type lists including
             * IQI, property, constitutive and stabilization dependencies
             * ( a list for master and a list for slave dof and dv types )
             */
            void build_global_dof_dv_and_field_type_list();

            //------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes  a cell of cell of dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dv_type_list(
                    const moris::Cell< moris::Cell< PDV_Type > >& aDvTypes,
                    mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< PDV_Type > >& get_dv_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * return a cell of field types active for the IQI
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] aFieldTypes a list of group of field types
             */
            const moris::Cell< moris::Cell< mtk::Field_Type > >& get_field_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

            //------------------------------------------------------------------------------
            /**
             * set IQI active field types
             * @param[ in ] aFieldTypes a list of group of field types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_field_type_list(
                    const moris::Cell< moris::Cell< mtk::Field_Type > >& aDvTypes,
                    mtk::Master_Slave                                    aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get a non unique list of dv type including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aGlobalDvTypeList a non unique list of dv types to fill
             */
            void get_non_unique_global_dv_type_list( moris::Cell< PDV_Type >& aGlobalDvTypeList );

            //------------------------------------------------------------------------------
            /**
             * get a unique global dv type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            const moris::Cell< moris::Cell< PDV_Type > >& get_global_dv_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get global field type list. TODO: Field types are only used by the IWG.
             * If a user wants to use them in a property or CM this cuntion has to be modified in the same way than get_global_dof_type_list()
             * @param[ in ]  aIsMaster    enum master or slave
             * @param[ out ] mFieldTypes global list of group of dv types
             */
            const moris::Cell< moris::Cell< mtk::Field_Type > >& get_global_field_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty       a property pointer
             * @param[ in ] aPropertyString a string describing the property
             * @param[ in ] aIsMaster       enum master or slave
             */
            void set_property(
                    std::shared_ptr< Property > aProperty,
                    std::string                 aPropertyString,
                    mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get master or slave properties
             * @param[ in ]  aIsMaster   enum master or slave
             * @param[ out ] aProperties cell of property pointers
             */
            moris::Cell< std::shared_ptr< fem::Property > >& get_properties(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set constitutive model
             * @param[ in ] aConstitutiveModel  a constitutive model pointer
             * @param[ in ] aConstitutiveString a string describing the constitutive model
             * @param[ in ] aIsMaster           an enum for master or slave
             */
            void set_constitutive_model(
                    std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                    std::string                           aConstitutiveString,
                    mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get master or slave constitutive models
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aConstitutiveModels cell of constitutive model pointers
             */
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > >& get_constitutive_models(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set stabilization parameter
             * @param[ in ] aStabilizationParameter a stabilization parameter pointer
             * @param[ in ] aStabilizationString    a string defining the stabilization parameter
             */
            void set_stabilization_parameter(
                    std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                    std::string                                aStabilizationString );

            //------------------------------------------------------------------------------
            /**
             * get stabilization parameters
             * @param[ out ] mStabilizationParam cell of stabilization parameter pointers
             */
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > >&
            get_stabilization_parameters()
            {
                // return stabilization parameter pointers
                return mStabilizationParam;
            }

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
                MORIS_ERROR( false, "IQI::compute_dQIdp - Not implemented for base class. " );
            }

            //------------------------------------------------------------------------------
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
                    moris::real                        aWStar,
                    moris::real                        aPerturbation,
                    fem::FDScheme_Type                 aFDSchemeType,
                    Matrix< DDSMat >&                  aGeoLocalAssembly,
                    moris::Cell< Matrix< IndexMat > >& aVertexIndices )
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
                    moris::real                        aWStar,
                    moris::real                        aPerturbation,
                    fem::FDScheme_Type                 aFDSchemeType,
                    Matrix< DDSMat >&                  aGeoLocalAssembly,
                    moris::Cell< Matrix< IndexMat > >& aVertexIndices );

            void select_dQIdp_FD_geometry_sideset(
                    moris::real                        aWStar,
                    moris::real                        aPerturbation,
                    fem::FDScheme_Type                 aFDSchemeType,
                    Matrix< DDSMat >&                  aGeoLocalAssembly,
                    moris::Cell< Matrix< IndexMat > >& aVertexIndices );

            void
            select_dQIdp_FD_geometry_double(
                    moris::real                        aWStar,
                    moris::real                        aPerturbation,
                    fem::FDScheme_Type                 aFDSchemeType,
                    Matrix< DDSMat >&                  aGeoLocalAssembly,
                    moris::Cell< Matrix< IndexMat > >& aVertexIndices )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdp_FD_geometry_double - not implemented yet" );
            }

            //------------------------------------------------------------------------------
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
             * Evaluate the quantity of interest.
             * @param[ in ] aWStar            weight associated to evaluation point
             */
            virtual void compute_QI( real aWStar ) = 0;

            //------------------------------------------------------------------------------
            /**
             * Compute the derivative of the quantities of interest wrt requested dof types.
             * @param[ in ]  aDofType Dof type being evaluated
             * @param[ out ] adQIdu derivative of quantity of interest
             */
            virtual void compute_dQIdu(
                    moris::Cell< MSI::Dof_Type >& aDofType,
                    Matrix< DDRMat >&             adQIdu ) = 0;

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
             * get the matrix dimension of the IQI in order to initialize the size of the mGloblaIQIVal
             * returns 1*1 which is a scaler by default
             */
            virtual std::pair< uint, uint >
            get_matrix_dim()
            {
                return std::make_pair( 1, 1 );
            }
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IQI_HPP_ */
