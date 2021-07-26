/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe/noel
 */
#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_
//MRS/CON/src
#include "cl_Cell.hpp"
//LNA/src
#include "cl_Matrix.hpp"
#include "typedefs.hpp"
#include "fn_vectorize.hpp"
#include "fn_isfinite.hpp"
//MRS/COR/src // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Enums.hpp"
#include "fn_FEM_FD_Scheme.hpp"
//FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
//GEN/src
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace fem
    {
        class Set;
        class Cluster;
        class Field_Interpolator_Manager;
        class FEM_Model;

        //------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
            protected :

                // FEM set pointer
                fem::Set * mSet = nullptr;

                // cluster pointer
                fem::Cluster * mCluster = nullptr;

                // nodal weak BCs
                Matrix< DDRMat > mNodalWeakBCs;

                // normal
                Matrix< DDRMat > mNormal;

                // residual dof type
                moris::Cell< moris::Cell< MSI::Dof_Type > > mResidualDofType;

                // master and slave dof type lists
                moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
                moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

                // bool for building global dof type list and map
                bool mGlobalDofBuild = true;
                bool mGlobalDvBuild = true;
                bool mGlobalFieldBuild = true;

                // master and slave global dof type lists
                moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
                moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

                // master and slave requested global dof type lists
                moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedMasterGlobalDofTypes;
                moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedSlaveGlobalDofTypes;

                // master and slave field interpolator managers
                Field_Interpolator_Manager * mMasterFIManager          = nullptr;
                Field_Interpolator_Manager * mSlaveFIManager           = nullptr;
                Field_Interpolator_Manager * mMasterPreviousFIManager  = nullptr;

                // master and slave dv type lists
                moris::Cell< moris::Cell< PDV_Type > > mMasterDvTypes;
                moris::Cell< moris::Cell< PDV_Type > > mSlaveDvTypes;

                // master and slave global dv type list
                moris::Cell< moris::Cell< PDV_Type > > mMasterGlobalDvTypes;
                moris::Cell< moris::Cell< PDV_Type > > mSlaveGlobalDvTypes;

                // master and slave field type lists
                moris::Cell< moris::Cell< mtk::Field_Type > > mMasterFieldTypes;
                moris::Cell< moris::Cell< mtk::Field_Type > > mSlaveFieldTypes;

                // master and slave global dv type list
                moris::Cell< moris::Cell< mtk::Field_Type > > mMasterGlobalFieldTypes;
                moris::Cell< moris::Cell< mtk::Field_Type > > mSlaveGlobalFieldTypes;

                // master and slave properties
                moris::Cell< std::shared_ptr< Property > > mMasterProp;
                moris::Cell< std::shared_ptr< Property > > mSlaveProp;

                // local string to int map for properties
                std::map< std::string, uint > mPropertyMap;

                // master and slave material models
                moris::Cell< std::shared_ptr< fem::Material_Model > > mMasterMM;
                moris::Cell< std::shared_ptr< fem::Material_Model > > mSlaveMM;

                // Local string to int map for material models
                std::map< std::string, uint > mMaterialMap;

                // master and slave constitutive models
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mMasterCM;
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mSlaveCM;

                // Local string to int map for constitutive models
                std::map< std::string, uint > mConstitutiveMap;

                // stabilization parameters
                moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > mStabilizationParam;

                // local string to int map for stabilizations
                std::map< std::string, uint > mStabilizationMap;

                // active cluster measure on IWG flag
                bool mActiveCMEAFlag = false;

                // interpolation order for IWG
                uint mOrder = MORIS_UINT_MAX;

                // bulk type
                fem::Element_Type mBulkType = fem::Element_Type::BULK;

                // strings for master and slave phase name
                std::string mMasterPhaseName;
                std::string mSlavePhaseName;

                // bool for time continuity
                bool mTimeContinuity = false;

                // bool for time boundary
                bool mTimeBoundary = false;

                // bool for ghost
                bool mIsGhost = false;

                // string for IWG name
                std::string mName;

                //! string for IWG name
                enum moris::fem::IWG_Type mIWGType = moris::fem::IWG_Type::UNDEFINED;

                // function pointers
                void ( IWG:: * m_compute_jacobian_FD )(
                        real               aWStar,
                        real               aPerturbation,
                        fem::FDScheme_Type aFDSchemeType,
                        bool               aUseAbsolutePerturbations ) = nullptr;
                void ( IWG:: * m_compute_dRdp_FD_material )(
                        real               aWStar,
                        real               aPerturbation,
                        fem::FDScheme_Type aFDSchemeType ) = nullptr;
                void ( IWG:: * m_compute_dRdp_FD_geometry )(
                        real                                aWStar,
                        real                                aPerturbation,
                        fem::FDScheme_Type                  aFDSchemeType,
                        Matrix< DDSMat >                  & aGeoLocalAssembly,
                        moris::Cell< Matrix< IndexMat > > & aVertexIndices ) = nullptr;

                // function pointer for building the perturbation size for FD
                real ( IWG:: * m_build_perturbation_size )(
                         const real & aPerturbation,
                         const real & aCoefficientToPerturb,
                         const real & aMaxPerturbation,
                         const real & aTolerance ) = nullptr;

                //------------------------------------------------------------------------------
            public :

                //------------------------------------------------------------------------------
                /**
                 * constructor
                 */
                IWG(){};

                //------------------------------------------------------------------------------
                /**
                 * virtual destructor
                 */
                virtual ~IWG(){};

                //------------------------------------------------------------------------------
                /**
                 * set name
                 * param[ in ] aName a string for CM name
                 */
                void set_name( std::string aName )
                {
                    mName = aName;
                }

                //------------------------------------------------------------------------------
                /**
                 * get name
                 * param[ out ] mName a string for CM name
                 */
                std::string get_name()
                {
                    return mName;
                }

                //------------------------------------------------------------------------------
                /**
                 * set time continuity flag
                 * param[ in ] aTimeContinuity bool true if IWG for time continuity
                 */
                void set_time_continuity( bool aTimeContinuity )
                {
                    mTimeContinuity = aTimeContinuity;
                }

                //------------------------------------------------------------------------------
                /**
                 * get time continuity flag
                 * param[ out ] mTimeContinuity ool true if IWG for time continuity
                 */
                bool get_time_continuity()
                {
                    return mTimeContinuity;
                }

                //------------------------------------------------------------------------------
                /**
                 * set time boundary flag
                 * param[ in ] aTimeBoundary bool true if IWG for time boundary
                 */
                void set_time_boundary( bool aTimeBoundary )
                {
                    mTimeBoundary = aTimeBoundary;
                }

                //------------------------------------------------------------------------------
                /**
                 * get time boundary flag
                 * param[ out ] mTimeBoundary bool true if IWG for time boundary
                 */
                bool get_time_boundary()
                {
                    return mTimeBoundary;
                }

                //------------------------------------------------------------------------------
                /**
                 * set ghost flag
                 * param[ in ] aIsGhost bool true if IWG for ghost
                 */
                void set_ghost_flag( bool aIsGhost )
                {
                    mIsGhost = aIsGhost;
                }

                //------------------------------------------------------------------------------
                /**
                 * get ghost flag
                 * param[ out ] mIsGhost bool true if IWG for ghost
                 */
                bool get_ghost_flag()
                {
                    return mIsGhost;
                }

                //------------------------------------------------------------------------------
                /**
                 * get IWG type
                 * param[ out ] mIWGType an enum of the IWG type. type only implemented for TIME_CONTINUITY_DOF. All others return UNDEFINED.
                 *              If needed you can implement the type for the others. Just folow the TIME_CONTINUITY_DOF
                 */
                enum moris::fem::IWG_Type get_IWG_type()
                {
                    return mIWGType;
                };

                //------------------------------------------------------------------------------
                /**
                 * set phase name
                 * param[ in ] aPhaseName a string for phase name
                 * param[ in ] aIsMaster  an enum for master or slave
                 */
                void set_phase_name(
                        std::string aPhaseName,
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * get phase name
                 * param[ in ]  aIsMaster an enum for master or slave
                 * param[ out ] mName     a string for phase name
                 */
                std::string get_phase_name( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * print name
                 */
                void print_names();

                //------------------------------------------------------------------------------
                /*
                 * set member set pointer
                 * @param[ in ] aSetPointer a FEM set pointer
                 */
                void set_set_pointer( Set * aSetPointer )
                {
                    mSet = aSetPointer;

                    // set function pointer
                    this->set_function_pointers();
                }

                //------------------------------------------------------------------------------
                /*
                 * set fem cluster pointer
                 * @param[ in ] aClusterPointer a FEM cluster pointer
                 */
                void set_cluster_pointer( fem::Cluster * aClusterPointer )
                {
                    mCluster = aClusterPointer;
                }

                //------------------------------------------------------------------------------
                /*
                 * set function pointers
                 */
                void set_function_pointers();

                //------------------------------------------------------------------------------
                /*
                 * get member set pointer
                 * @param[ out ] aSetPointer a FEM set pointer
                 */
                Set * get_set_pointer()
                {
                    return mSet;
                }

                //------------------------------------------------------------------------------
                /*
                 * set field interpolator manager
                 * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
                 * @param[ in ] aIsMaster                 an enum for master or slave
                 */
                void set_field_interpolator_manager(
                        Field_Interpolator_Manager * aFieldInterpolatorManager,
                        mtk::Master_Slave            aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /*
                 * set field interpolator manager for previous time step
                 * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
                 * @param[ in ] aIsMaster                 an enum for master or slave
                 */
                void set_field_interpolator_manager_previous_time(
                        Field_Interpolator_Manager * aFieldInterpolatorManager,
                        mtk::Master_Slave            aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /*
                 * get field interpolator manager
                 * @param[ out ] aFieldInterpolatorManager a field interpolator manager pointer
                 * @param[ in ]  aIsMaster                 an enum for master or slave
                 */
                Field_Interpolator_Manager * get_field_interpolator_manager(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /*
                 * free memory
                 */
                void free_memory(){}

                //------------------------------------------------------------------------------
                /**
                 * set nodal weak BCs
                 * @param[ in ] aNodalWeakBCs matrix with nodal values
                 */
                void set_nodal_weak_bcs( Matrix< DDRMat > & aNodalWeakBCs )
                {
                    mNodalWeakBCs = aNodalWeakBCs;
                }

                //------------------------------------------------------------------------------
                /**
                 * set normal
                 * @param[ in ] aNormal normal vector
                 */
                void set_normal( Matrix< DDRMat > & aNormal );

                //------------------------------------------------------------------------------
                /**
                 * set residual dof type
                 * @param[ in ] aResidualdofType a cell of residual dof types
                 */
                void set_residual_dof_type( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aResidualDofType )
                {
                    mResidualDofType = aResidualDofType;
                }

                //------------------------------------------------------------------------------
                /**
                 * return a dof type for the residual
                 * @param[ out ] aResidualdofType a cell of residual dof types
                 */
                const moris::Cell< moris::Cell < MSI::Dof_Type > > & get_residual_dof_type() const
                {
                    return mResidualDofType;
                }

                //------------------------------------------------------------------------------
                /**
                 * set interpolation order for the residual dof type
                 */
                void set_interpolation_order();

                void set_interpolation_order( uint aOrder )
                {
                    // set order
                    mOrder = aOrder;
                }

                //------------------------------------------------------------------------------
                /**
                 * set bulk type
                 * @param[ in ] aBulkType element type for the IWG
                 */
                void set_bulk_type( fem::Element_Type aBulkType )
                {
                    mBulkType = aBulkType;
                }

                //------------------------------------------------------------------------------
                /**
                 * get bulk type
                 * @param[ out ] mBulkType element type for the IWG
                 */
                fem::Element_Type get_bulk_type()
                {
                    return mBulkType;
                }

                //------------------------------------------------------------------------------
                /**
                 * set IWG active dof types
                 * @param[ in ] aDofTypes a list of group of dof types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dof_type_list(
                        const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        mtk::Master_Slave                                   aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * return a cell of dof types active for the IWG
                 * @param[ in ] aIsMaster enum master or slave
                 * @param[ out ] aDofTypes a list of group of dof types
                 */
                const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //------------------------------------------------------------------------------
                /**
                 * set IWG active dv types
                 * @param[ in ] aDvTypes a list of group of dv types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_dv_type_list(
                        const moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                              mtk::Master_Slave                        aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * return a cell of dv types active for the IWG
                 * @param[ in ]  aIsMaster enum master or slave
                 * @param[ out ] aDvTypes a list of group of dv types
                 */
                const moris::Cell< moris::Cell< PDV_Type > > & get_dv_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //------------------------------------------------------------------------------
                /**
                 * return a cell of field types active for the IWG
                 * @param[ in ]  aIsMaster enum master or slave
                 * @param[ out ] aFieldTypes a list of group of field types
                 */
                const moris::Cell< moris::Cell< mtk::Field_Type > > & get_field_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const;

                //------------------------------------------------------------------------------
                /**
                 * set IWG active field types
                 * @param[ in ] aFieldTypes a list of group of field types
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                void set_field_type_list(
                        const moris::Cell< moris::Cell< mtk::Field_Type > > & aDvTypes,
                              mtk::Master_Slave                        aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * check that field interpolators were assigned
                 * @param[ in ]  aIsMaster enum master or slave
                 */
                void check_field_interpolators(
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
                 * get properties
                 * @param[ in ]  aIsMaster   enum master or slave
                 * @param[ out ] aProperties cell of property pointers
                 */
                moris::Cell< std::shared_ptr< Property > > & get_properties(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * set material model
                 * @param[ in ] aMaterialModel       a material model pointer
                 * @param[ in ] aMaterialModelString a string defining the material model
                 * @param[ in ] aIsMaster            an enum for master or slave
                 */
                void set_material_model(
                        std::shared_ptr< Material_Model > aMaterialModel,
                        std::string                       aMaterialModelString,
                        mtk::Master_Slave                 aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * get material models
                 * @param[ in ]  aIsMaster           enum master or slave
                 * @param[ out ] aMaterialModels     cell of material model pointers
                 */
                moris::Cell< std::shared_ptr< Material_Model > > & get_material_models(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * set constitutive model
                 * @param[ in ] aConstitutiveModel  a constitutive model pointer
                 * @param[ in ] aConstitutiveString a string defining the constitutive model
                 * @param[ in ] aIsMaster           an enum for master or slave
                 */
                void set_constitutive_model(
                        std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                        std::string                           aConstitutiveString,
                        mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * get constitutive models
                 * @param[ in ]  aIsMaster           enum master or slave
                 * @param[ out ] aConstitutiveModels cell of constitutive model pointers
                 */
                moris::Cell< std::shared_ptr< Constitutive_Model > > & get_constitutive_models(
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
                moris::Cell< std::shared_ptr< Stabilization_Parameter > > & get_stabilization_parameters()
                {
                    // return penalty parameter pointers
                    return mStabilizationParam;
                }

                //------------------------------------------------------------------------------
                /**
                 * create a global dof type list including
                 * IWG, property, constitutive and stabilization dependencies
                 */
                void build_global_dof_dv_and_field_type_list();

                //------------------------------------------------------------------------------
                /**
                 * get a non unique list of dof type including
                 * IWG, property, constitutive and stabilization dependencies
                 * for both master and slave
                 */
                void get_non_unique_dof_dv_and_field_types(
                        moris::Cell< moris::Cell< MSI::Dof_Type > >   & aDofTypes,
                        moris::Cell< moris::Cell< PDV_Type > >        & aDvTypes,
                        moris::Cell< moris::Cell< mtk::Field_Type > > & aFieldTypes );

                //------------------------------------------------------------------------------
                /**
                 * get global dof type list
                 * @param[ in ]  aIsMaster       enum master or slave
                 * @param[ out ] mGlobalDofTypes global list of group of dof types
                 */
                const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * get global dv type list
                 * @param[ in ]  aIsMaster       enum master or slave
                 * @param[ out ] mGlobalDvTypes global list of group of dv types
                 */
                const moris::Cell< moris::Cell< PDV_Type > > & get_global_dv_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * get global field type list. TODO: Field types are only used by the IWG.
                 * If a user wants to use them in a property or CM this cuntion has to be modified in the same way than get_global_dof_type_list()
                 * @param[ in ]  aIsMaster    enum master or slave
                 * @param[ out ] mFieldTypes global list of group of dv types
                 */
                const moris::Cell< moris::Cell< mtk::Field_Type > > & get_global_field_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the residual
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_residual( real aWStar ) = 0;

                //------------------------------------------------------------------------------
                /**
                 * evaluate the Jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_jacobian( real aWStar ) = 0;

                virtual void compute_jacobian_previous( real aWStar )
                {
                    MORIS_ERROR(false, "compute_jacobian_previous() not implemented");
                };

                //------------------------------------------------------------------------------
                /**
                 * evaluate the Jacobian by finite difference
                 * @param[ in ] aPerturbation real to perturb for FD
                 * @param[ in ] aWStar        weight associated to evaluation point
                 * @param[ in ] aFDSchemeType enum for FD scheme
                 */
                void compute_jacobian_FD(
                        real               aWStar,
                        real               aPerturbation,
                        fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_5,
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

                //------------------------------------------------------------------------------
                /**
                 * evaluate the residual and the Jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                virtual void compute_jacobian_and_residual( real aWStar ) = 0;

                //------------------------------------------------------------------------------
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
                        real               aPerturbation,
                        real               aEpsilon,
                        real               aWStar,
                        Matrix< DDRMat > & aJacobians,
                        Matrix< DDRMat > & aJacobiansFD,
                        bool               aErrorPrint = false,
                        bool               aUseAbsolutePerturbations = false );

                //------------------------------------------------------------------------------
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
                        real               aPerturbation,
                        real               aEpsilon,
                        real               aWStar,
                        Matrix< DDRMat > & aJacobians,
                        Matrix< DDRMat > & aJacobiansFD,
                        bool               aErrorPrint = false,  
                        bool               aMaxErrorPrint = false,
                        moris::real        aFDtolerance = -1.0,
                        bool               aUseAbsolutePerturbations = false );                        

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the residual wrt the design variables
                 * @param[ in ] aWStar weight associated to evaluation point
                 */
                virtual void compute_dRdp( real aWStar ) = 0;

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the residual
                 * wrt the material design variables by finite difference
                 * @param[ in ] aWStar        weight associated to evaluation point
                 * @param[ in ] aPerturbation real for dv perturbation
                 * @param[ in ] aFDSchemeType enum for FD scheme
                 */
                void compute_dRdp_FD_material(
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

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the residual
                 * wrt the geometry design variables by finite difference
                 * @param[ in ] aWStar            weight associated to evaluation point
                 * @param[ in ] aPerturbation     real for relative dv perturbation
                 * @param[ in ] aGeoLocalAssembly matrix filled with pdv local assembly indices
                 * @param[ in ] aFDSchemeType     enum for FD scheme
                 */
                void compute_dRdp_FD_geometry(
                        moris::real                         aWStar,
                        moris::real                         aPerturbation,
                        fem::FDScheme_Type                  aFDSchemeType,
                        Matrix< DDSMat >                  & aGeoLocalAssembly,
                        moris::Cell< Matrix< IndexMat > > & aVertexIndices )
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
                        moris::real                         aWStar,
                        moris::real                         aPerturbation,
                        fem::FDScheme_Type                  aFDSchemeType,
                        Matrix< DDSMat >                  & aGeoLocalAssembly,
                        moris::Cell< Matrix< IndexMat > > & aVertexIndices );

                void select_dRdp_FD_geometry_sideset(
                        moris::real                         aWStar,
                        moris::real                         aPerturbation,
                        fem::FDScheme_Type                  aFDSchemeType,
                        Matrix< DDSMat >                  & aGeoLocalAssembly,
                        moris::Cell< Matrix< IndexMat > > & aVertexIndices );

                void select_dRdp_FD_geometry_time_sideset(
                        moris::real                         aWStar,
                        moris::real                         aPerturbation,
                        fem::FDScheme_Type                  aFDSchemeType,
                        Matrix< DDSMat >                  & aGeoLocalAssembly,
                        moris::Cell< Matrix< IndexMat > > & aVertexIndices );

                void select_dRdp_FD_geometry_double(
                        moris::real                         aWStar,
                        moris::real                         aPerturbation,
                        fem::FDScheme_Type                  aFDSchemeType,
                        Matrix< DDSMat >                  & aGeoLocalAssembly,
                        moris::Cell< Matrix< IndexMat > > & aVertexIndices );

                //------------------------------------------------------------------------------
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

                //------------------------------------------------------------------------------
                /**
                 * build perturbation size for finite difference
                 * @param[ in ] aPerturbation         provided perturbation size from input
                 * @param[ in ] aCoefficientToPerturb coefficient to perturb
                 * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
                 */
                real build_perturbation_size(
                        const real & aPerturbation,
                        const real & aCoefficientToPerturb,
                        const real & aMaxPerturbation,
                        const real   aTolerance = 1e-12 );

                /**
                 * build relative perturbation size for finite difference
                 * @param[ in ] aPerturbation         provided perturbation size from input
                 * @param[ in ] aCoefficientToPerturb coefficient to perturb
                 * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
                 */
                real build_perturbation_size_relative(
                        const real & aPerturbation,
                        const real & aCoefficientToPerturb,
                        const real & aMaxPerturbation,
                        const real & aTolerance );

                /**
                 * build absolute perturbation size for finite difference
                 * @param[ in ] aPerturbation         provided perturbation size from input
                 * @param[ in ] aCoefficientToPerturb coefficient to perturb
                 * @param[ in ] aTolerance            tolerance to check that built perturbation is not too small
                 */
                real build_perturbation_size_absolute(
                        const real & aPerturbation,
                        const real & aCoefficientToPerturb,
                        const real & aMaxPerturbation,
                        const real & aTolerance );

                //------------------------------------------------------------------------------
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
                        		const real & aPerturbation,
                                const real & aCoefficientToPerturb,
								const uint & aSpatialDirection,
                                fem::FDScheme_Type & aUsedFDScheme );

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags specific to child IWG
                 */
                virtual void reset_spec_eval_flags(){};

                //------------------------------------------------------------------------------
                /**
                 * build a list of dof types requested by the solver and owned by the IWG
                 * @param[ in ] aIsResidual bool true if residual evaluation
                 */
                void build_requested_dof_type_list( const bool aIsStaggered );

        };
        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
