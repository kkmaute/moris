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
#include "fn_reshape.hpp"
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
        class Field_Interpolator_Manager;
        //------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
            protected :

                // FEM set pointer
                fem::Set * mSet = nullptr;

                // nodal weak BCs
                Matrix< DDRMat > mNodalWeakBCs;

                // normal
                Matrix< DDRMat > mNormal;

                // residual dof type
                moris::Cell< MSI::Dof_Type > mResidualDofType;

                // master and slave dof type lists
                moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
                moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

                // bool for building global dof type list and map
                bool mGlobalDofBuild = true;
                bool mGlobalDvBuild = true;

                // master and slave global dof type lists
                moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
                moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

                // master and slave requested global dof type lists
                moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedMasterGlobalDofTypes;
                moris::Cell< moris::Cell< MSI::Dof_Type > > mRequestedSlaveGlobalDofTypes;

                // master and slave field interpolator managers
                Field_Interpolator_Manager * mMasterFIManager = nullptr;
                Field_Interpolator_Manager * mSlaveFIManager  = nullptr;
                Field_Interpolator_Manager * mMasterPreviousFIManager  = nullptr;

                // master and slave dv type lists
                moris::Cell< moris::Cell< PDV_Type > > mMasterDvTypes;
                moris::Cell< moris::Cell< PDV_Type > > mSlaveDvTypes;

                // master and slave global dv type list
                moris::Cell< moris::Cell< PDV_Type > > mMasterGlobalDvTypes;
                moris::Cell< moris::Cell< PDV_Type > > mSlaveGlobalDvTypes;

                // master and slave properties
                moris::Cell< std::shared_ptr< Property > > mMasterProp;
                moris::Cell< std::shared_ptr< Property > > mSlaveProp;

                // master and slave constitutive models
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mMasterCM;
                moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mSlaveCM;

                // stabilization parameters
                moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > mStabilizationParam;

                // string for IWG name
                std::string mName;

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
                }

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
                void set_residual_dof_type( const moris::Cell< MSI::Dof_Type > & aResidualDofType )
                {
                    mResidualDofType = aResidualDofType;
                }

                //------------------------------------------------------------------------------
                /**
                 * return a dof type for the residual
                 * @param[ out ] aResidualdofType a cell of residual dof types
                 */
                const moris::Cell< MSI::Dof_Type > & get_residual_dof_type() const
                {
                    return mResidualDofType;
                };

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
                virtual void set_property(
                        std::shared_ptr< Property > aProperty,
                        std::string                 aPropertyString,
                        mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER )
                {
                    MORIS_ASSERT( false, "IWG::set_property - This function does nothing.");
                }

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
                 * set constitutive model
                 * @param[ in ] aConstitutiveModel  a constitutive model pointer
                 * @param[ in ] aConstitutiveString a string defining the constitutive model
                 * @param[ in ] aIsMaster           an enum for master or slave
                 */
                virtual void set_constitutive_model(
                        std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                        std::string                           aConstitutiveString,
                        mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
                {
                    MORIS_ERROR( false, "IWG::set_constitutive_model - This function does nothing." );
                }

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
                virtual void set_stabilization_parameter(
                        std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                        std::string                                aStabilizationString )
                {
                    MORIS_ERROR( false, "IWG::set_stabilization_parameter - This function does nothing." );
                }

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
                void build_global_dof_and_dv_type_list();

                //------------------------------------------------------------------------------
                /**
                 * get a non unique list of dof type including
                 * IWG, property, constitutive and stabilization dependencies
                 * for both master and slave
                 */
                void get_non_unique_dof_and_dv_types(
                        moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                        moris::Cell< moris::Cell< PDV_Type > >      & aDvTypes );

                //------------------------------------------------------------------------------
                /**
                 * get global dof type list
                 * @param[ in ]  aIsMaster       enum master or slave
                 * @param[ out ] mGlobalDofTypes global list of group of dof types
                 */
                moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list(
                        mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
                /**
                 * get global dv type list
                 * @param[ in ]  aIsMaster       enum master or slave
                 * @param[ out ] mGlobalDvTypes global list of group of dv types
                 */
                moris::Cell< moris::Cell< PDV_Type > > & get_global_dv_type_list(
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
                        fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_5 );

                void compute_jacobian_FD_double(
                        real               aWStar,
                        real               aPerturbation,
                        fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_5 );

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
                        bool               aErrorPrint = false );

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
                 */
                void compute_dRdp_FD_material(
                        moris::real        aWStar,
                        moris::real        aPerturbation,
                        fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

                void compute_dRdp_FD_material_double(
                        moris::real        aWStar,
                        moris::real        aPerturbation,
                        fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

                //------------------------------------------------------------------------------
                /**
                 * evaluate the derivative of the residual
                 * wrt the geometry design variables by finite difference
                 * @param[ in ] aWStar         weight associated to evaluation point
                 * @param[ in ] aPerturbation  real for dv perturbation
                 * @param[ in ] aIsActive      cell of vectors for active dv
                 * @param[ in ] aVertexIndices vertices indices
                 */
                void compute_dRdp_FD_geometry(
                        moris::real                       aWStar,
                        moris::real                       aPerturbation,
                        moris::Cell< Matrix< DDSMat > > & aIsActive,
                        Matrix< IndexMat >              & aVertexIndices,
                        fem::FDScheme_Type                aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

                void compute_dRdp_FD_geometry_double(
                        moris::real                       aWStar,
                        moris::real                       aPerturbation,
                        moris::Cell< Matrix< DDSMat > > & aMasterIsActive,
                        Matrix< IndexMat >              & aMasterVertexIndices,
                        moris::Cell< Matrix< DDSMat > > & aSlaveIsActive,
                        Matrix< IndexMat >              & aSlaveVertexIndices,
                        fem::FDScheme_Type                aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

                //------------------------------------------------------------------------------
                /**
                 * reset evaluation flags
                 */
                void reset_eval_flags();

                //------------------------------------------------------------------------------
                /**
                 * build a list of dof types requested by the solver and owned by the IWG
                 * @param[ in ] aIsResidual bool true if residual evaluation
                 */
                void build_requested_dof_type_list( const bool aIsResidual );

        };
        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
