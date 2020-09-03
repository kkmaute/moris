/*
 * cl_FEM_IQI.hpp
 *
 *  Created on: Oct 08, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_IQI_HPP_
#define SRC_FEM_CL_FEM_IQI_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
//MRS/COR/src           // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
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
//FEM/VIS/src
#include "cl_VIS_Output_Enums.hpp"
//GEN/src
#include "cl_GEN_Pdv_Enums.hpp"
//LINALG/src
#include "fn_reshape.hpp"

namespace moris
{
    namespace fem
    {
        class Set;
        class Field_Interpolator_Manager;
        //------------------------------------------------------------------------------
        /**
         * Integrand of a quantity of interest
         */
        class IQI
        {
        protected :

            // FEM set pointer
            fem::Set * mSet = nullptr;

            // IQI type
            enum vis::Output_Type mIQIType;

            // FEM IQI type
            enum fem::IQI_Type mFEMIQIType;

            // Phase type
            enum Phase_Type mIQIMatType;

            // IQI type index
            sint mIQITypeIndex = -1;

            // normal
            Matrix< DDRMat > mNormal;

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
            Field_Interpolator_Manager * mMasterFIManager = nullptr;
            Field_Interpolator_Manager * mSlaveFIManager  = nullptr;

            // master and slave dv type lists
            moris::Cell< moris::Cell< PDV_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< PDV_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< PDV_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< PDV_Type > > mSlaveGlobalDvTypes;

            // flag for building global dv type list
            bool mGlobalDvBuild = true;

            // master and slave properties
            moris::Cell< std::shared_ptr< fem::Property > > mMasterProp;
            moris::Cell< std::shared_ptr< fem::Property > > mSlaveProp;

            // master and slave constitutive models
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mMasterCM;
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > mSlaveCM;

            // stabilization parameters
            moris::Cell< std::shared_ptr< Stabilization_Parameter > > mStabilizationParam;

            // local string to dof enum map
            std::map< std::string, MSI::Dof_Type > mMasterDofMap;
            std::map< std::string, MSI::Dof_Type > mSlaveDofMap;

            // parameters
            moris::Cell< Matrix< DDRMat > > mParameters;

            // IQI name
            std::string mName;

        private:

            // Normalization
            real mReferenceValue = 1.0;
            bool mNormalized = false;

            //------------------------------------------------------------------------------
        public :

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
             * print names
             */
            void print_names();

            //------------------------------------------------------------------------------
            /**
             * get vis IQI type
             */
            enum vis::Output_Type get_IQI_type()
            {
                return mIQIType;
            }

            //------------------------------------------------------------------------------
            /**
             * get fem IQI type
             */
            enum fem::IQI_Type get_fem_IQI_type()
            {
                return mFEMIQIType;
            }

            //------------------------------------------------------------------------------
            /**
             * get IQI mat type
             */
            enum Phase_Type get_IQI_phase_type()
            {
                return mIQIMatType;
            }

            /**
             * set IQI mat type
             */
            void set_IQI_phase_type( enum Phase_Type aMatType )
            {
                mIQIMatType = aMatType;
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
            void set_reference_value(real aReferenceValue);

            //------------------------------------------------------------------------------
            /*
             * set fem set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_set_pointer( Set * aSetPointer )
            {
                mSet = aSetPointer;
            }

            //------------------------------------------------------------------------------
            /*
             * set output type
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_output_type( enum vis::Output_Type aOutputType )
            {
                mIQIType = aOutputType;
            }

            //------------------------------------------------------------------------------
            /*
             * set output type index
             * @param[ in ] aOutputTypeIndex output type index
             */
            void set_output_type_index( sint aOutputTypeIndex )
            {
                mIQITypeIndex = aOutputTypeIndex;
            }

            //------------------------------------------------------------------------------
            /**
             * set normal
             * @param[ in ] aNormal normal vector
             */
            void set_normal( Matrix< DDRMat > & aNormal )
            {
                mNormal = aNormal;
            }

            //------------------------------------------------------------------------------
            /**
             * set parameters
             * @param[ in ] aParameters a list of parameters
             */
            void set_parameters( moris::Cell< Matrix< DDRMat > > aParameters )
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
            /**
             * set dof type list for IQI
             * @param[ in ] aDofTypes a cell of group of dof types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dof_type_list(
                    const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                    mtk::Master_Slave                                   aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes   list of group of dof types
             * @param[ in ] aDofStrings list of strings describing the group dof types
             * @param[ in ] aIsMaster   enum for master or slave
             */
            virtual void set_dof_type_list(
                    moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                    moris::Cell< std::string >                  & aDofStrings,
                    mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER )
            {
                MORIS_ERROR( false, "IQI::set_dof_type_list - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get a non unique list of dof type including
             * IQI, property, constitutive and stabilization dependencies
             * for both master and slave
             */
            void get_non_unique_dof_and_dv_types(
                    moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                    moris::Cell< moris::Cell< PDV_Type > >      & aDvTypes );

            //------------------------------------------------------------------------------
            /**
             * get a unique global dof type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * build a global dof and dv type lists including
             * IQI, property, constitutive and stabilization dependencies
             * ( a list for master and a list for slave dof and dv types )
             */
            void build_global_dof_and_dv_type_lists();

            //------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes  a cell of cell of dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dv_type_list(
                    const moris::Cell< moris::Cell< PDV_Type > > & aDvTypes,
                    mtk::Master_Slave                              aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< PDV_Type > > & get_dv_type_list(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * get a non unique list of dv type including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aGlobalDvTypeList a non unique list of dv types to fill
             */
            void get_non_unique_global_dv_type_list( moris::Cell< PDV_Type > & aGlobalDvTypeList );

            //------------------------------------------------------------------------------
            /**
             * get a unique global dv type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< PDV_Type > > & get_global_dv_type_list(
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
                MORIS_ERROR( false, "IQI::set_property - This function does nothing.");
            }

            //------------------------------------------------------------------------------
            /**
             * get master or slave properties
             * @param[ in ]  aIsMaster   enum master or slave
             * @param[ out ] aProperties cell of property pointers
             */
            moris::Cell< std::shared_ptr< fem::Property > > & get_properties(
                    mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            //------------------------------------------------------------------------------
            /**
             * set constitutive model
             * @param[ in ] aConstitutiveModel  a constitutive model pointer
             * @param[ in ] aConstitutiveString a string describing the constitutive model
             * @param[ in ] aIsMaster           an enum for master or slave
             */
            virtual void set_constitutive_model(
                    std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                    std::string                           aConstitutiveString,
                    mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
            {
                MORIS_ERROR( false, "IQI::set_constitutive_model - This function does nothing." );
            }

            //------------------------------------------------------------------------------
            /**
             * get master or slave constitutive models
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aConstitutiveModels cell of constitutive model pointers
             */
            moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & get_constitutive_models(
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
                MORIS_ERROR( false, "IQI::set_stabilization_parameter - This function does nothing." );
            }

            //------------------------------------------------------------------------------
            /**
             * get stabilization parameters
             * @param[ out ] mStabilizationParam cell of stabilization parameter pointers
             */
            moris::Cell< std::shared_ptr< fem::Stabilization_Parameter > > & get_stabilization_parameters()
            {
                // return stabilization parameter pointers
                return mStabilizationParam;
            }

            //------------------------------------------------------------------------------

            /**
             * Gets the quantity of interest, including scaling.
             * @param[ in ] aQIVal quantity of interest matrix to fill
             */
            void get_QI( Matrix< DDRMat > & aQIVal );

            /**
             * Gets the derivative of the quantities of interest wrt requested dof types, including scaling.
             * @param[ in ] aDofType Dof type being evaluated
             * @param[ out ] adQIdu derivative of quantity of interest
             */
            void get_dQIdu( MSI::Dof_Type aDofType, Matrix< DDRMat > & adQIdu );

            //------------------------------------------------------------------------------
            /**
             * compute the quantities of interest, and add to the cell contained on the set
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void add_QI_on_set( real aWStar );

            /**
             * compute the derivatives of the quantities of interest, and add to the cell contained on the set
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void add_dQIdu_on_set( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * get requested dof type
             * @param[ in ] mRequestedDofType list of requested dof type
             */
            moris::Cell < enum MSI::Dof_Type > get_requested_dof_types();

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
            void compute_dQIdu_FD(
                    real               aWStar,
                    real               aPerturbation,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

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
                    Matrix< DDRMat > & adQIdu,
                    Matrix< DDRMat > & adQIduFD,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest wrt to dv types
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_dQIdp( real aWStar )
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
            void compute_dQIdp_FD_material(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

            void compute_dQIdp_FD_material_double(
                    moris::real        aWStar,
                    moris::real        aPerturbation,
                    fem::FDScheme_Type aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdp_FD_material_double - not implemented yet");
            }

            //------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest
             * wrt to geometry dv by finite difference
             * @param[ in ] aWStar         weight associated to evaluation point
             * @param[ in ] aPerturbation  dv relative perturbation
             * @param[ in ] aIsActive      cell of vectors for active dv
             * @param[ in ] aVertexIndices vertices indices
             * @param[ in ] aFDSchemeType enum for FD scheme
             */
            void compute_dQIdp_FD_geometry(
                    moris::real                       aWStar,
                    moris::real                       aPerturbation,
                    moris::Cell< Matrix< DDSMat > > & aIsActive,
                    Matrix< IndexMat >              & aVertexIndices,
                    fem::FDScheme_Type                aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL );

            void compute_dQIdp_FD_geometry_double(
                    moris::real                       aWStar,
                    moris::real                       aPerturbation,
                    moris::Cell< Matrix< DDSMat > > & aMasterIsActive,
                    Matrix< IndexMat >              & aMasterVertexIndices,
                    moris::Cell< Matrix< DDSMat > > & aSlaveIsActive,
                    Matrix< IndexMat >              & aSlaveVertexIndices,
                    fem::FDScheme_Type                aFDSchemeType = fem::FDScheme_Type::POINT_3_CENTRAL )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdp_FD_geometry_double - not implemented yet");
            }

        private:

            /**
             * Evaluate the quantity of interest.
             * @param[ out ] aQIVal quantity of interest matrix to fill
             */
            virtual void compute_QI( Matrix< DDRMat > & aQIVal ) = 0;

            //------------------------------------------------------------------------------
            /**
             * Compute the derivative of the quantities of interest wrt requested dof types.
             * @param[ in ] aDofType Dof type being evaluated
             * @param[ out ] adQIdu derivative of quantity of interest
             */
            virtual void compute_dQIdu( MSI::Dof_Type aDofType, Matrix< DDRMat > & adQIdu ) = 0;

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IQI_HPP_ */
