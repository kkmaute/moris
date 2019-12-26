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
#include "linalg_typedefs.hpp"              //MRS/COR/src           // note: linalg_typedefs.hpp must be included AFTER the cl_Matrix.hpp
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Stabilization_Parameter.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/MSI/src
#include "cl_VIS_Output_Enums.hpp"

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

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // flag for building global dof type list
            bool mGlobalDofBuild = true;

            // field interpolator manager pointer
            Field_Interpolator_Manager * mMasterFIManager = nullptr;
            Field_Interpolator_Manager * mSlaveFIManager  = nullptr;

            // master and slave dv type lists
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveGlobalDvTypes;

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
             * get IQI type
             */
            enum vis::Output_Type get_IQI_type()
            {
                return mIQIType;
            }

//------------------------------------------------------------------------------
            /**
             * rest evaluation flags for the IQI
             */
            void reset_eval_flags();

//------------------------------------------------------------------------------
            /*
             * set fel set pointer
             * @param[ in ] aSetPointer a FEM set pointer
             */
            void set_set_pointer( Set * aSetPointer )
            {
                mSet = aSetPointer;
            }

//------------------------------------------------------------------------------
            /*
             * set field interpolator manager pointer
             * @param[ in ] aFieldInterpolatorManager a field interpolator manager pointer
             */
            void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
                                                 mtk::Master_Slave            aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /*
             * get field interpolator manager
             * @param[ out ] aFieldInterpolatorManager a field interpolator manager pointer
             * @param[ in ]  aIsMaster                 an enum for master or slave
             */
            Field_Interpolator_Manager * get_field_interpolator_manager( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ) :
                        return mMasterFIManager;

                    case ( mtk::Master_Slave::SLAVE ) :
                        return mSlaveFIManager;

                    default :
                    {
                        MORIS_ERROR( false, "IWG::get_field_inetrpolator_manager - can only be master or slave." );
                        return mMasterFIManager;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dof_type_list( const moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                                    mtk::Master_Slave                                   aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterDofTypes = aDofTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlaveDofTypes = aDofTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "IQI::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterDofTypes;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveDofTypes;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IQI::get_dof_type_list - can only be master or slave." );
                        return mMasterDofTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * get a non unique list of dof type including
             * IQI, property, constitutive and stabilization dependencies
             * for both master and slave
             */
            void get_non_unique_global_dof_type_list( moris::Cell< MSI::Dof_Type > & aDofTypes );

//------------------------------------------------------------------------------
            /**
             * get a unique global dof type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // if the global list was not yet built
                if( mGlobalDofBuild )
                {
                    // build global dof type list
                    this->build_global_dof_type_list();

                    // update build flag
                    mGlobalDofBuild = false;
                }

                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterGlobalDofTypes;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveGlobalDofTypes;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IQI::get_global_dof_type_list - can only be master or slave." );
                        return mMasterGlobalDofTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * build a global dof type list including
             * IQI, property, constitutive and stabilization dependencies
             * ( a list for master and a list for slave do types )
             */
            void build_global_dof_type_list();

////------------------------------------------------------------------------------
//            /**
//             * set master or slave dof field interpolators for the IQI
//             * properties, constitutive models and stabilization parameters
//             * @param[ in ] aIsMaster an enum for master or slave
//             */
//            void set_dof_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes  a cell of cell of dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dv_type_list( const moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes,
                                   mtk::Master_Slave                                  aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterDvTypes = aDvTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlaveDvTypes = aDvTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "IQI::set_dv_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dv type list
                        return mMasterDvTypes;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dv type list
                        return mSlaveDvTypes;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IQI::get_dv_type_list - can only be master or slave." );
                        return mMasterDvTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * get a non unique list of dv type including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aGlobalDvTypeList a non unique list of dv types to fill
             */
            void get_non_unique_global_dv_type_list( moris::Cell< MSI::Dv_Type > & aGlobalDvTypeList );

//------------------------------------------------------------------------------
            /**
             * get a unique global dv type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // if the global list was not yet built
                if( mGlobalDvBuild )
                {
                    // build global dv type list
                    this->build_global_dv_type_list();

                    // update build flag
                    mGlobalDvBuild = false;
                }

                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dv type list
                        return mMasterGlobalDvTypes;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dv type list
                        return mSlaveGlobalDvTypes;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IQI::get_global_dv_type_list - can only be master or slave." );
                        return mMasterGlobalDvTypes;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * build a global dv type list including
             * IQI, property, constitutive and stabilization dependencies
             */
            void build_global_dv_type_list();

//------------------------------------------------------------------------------
            /**
             * set master or slave dv field interpolators for the IQI
             * properties, constitutive models and stabilization parameters
             * @param[ in ] aIsMaster an enum for master or slave
             */
            void set_dv_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /**
             * set geometry interpolator
             * @param[ in ] aGeometryInterpolator geometry interpolator pointers
             * @param[ in ] aIsMaster             enum for master or slave
             */
            void set_geometry_interpolator( Geometry_Interpolator* aGeometryInterpolator,
                                            mtk::Master_Slave      aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /**
             * set property
             * @param[ in ] aProperty       a property pointer
             * @param[ in ] aPropertyString a string describing the property
             * @param[ in ] aIsMaster       enum master or slave
             */
            virtual void set_property( std::shared_ptr< Property > aProperty,
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
             moris::Cell< std::shared_ptr< fem::Property > > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterProp;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveProp;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "IQI::get_properties - can only be master or slave." );
                         return mMasterProp;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * set constitutive model
              * @param[ in ] aConstitutiveModel  a constitutive model pointer
              * @param[ in ] aConstitutiveString a string describing the constitutive model
              * @param[ in ] aIsMaster           an enum for master or slave
              */
             virtual void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
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
             moris::Cell< std::shared_ptr< fem::Constitutive_Model > > & get_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterCM;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveCM;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "IQI::get_constitutive_models - can only be master or slave." );
                         return mMasterCM;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * set stabilization parameter
              * @param[ in ] aStabilizationParameter a stabilization parameter pointer
              * @param[ in ] aStabilizationString    a string defining the stabilization parameter
              */
             virtual void set_stabilization_parameter( std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
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
             * evaluate the quantity of interest
             * @param[ in ] aQIVal quantity of interest matrix to fill
             */
            virtual void compute_QI( Matrix< DDRMat > & aQIVal )
            {
                MORIS_ERROR( false, "IQI::compute_QI - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest wrt to dof types
             * @param[ in ] adQIdDof matrix to fill with derivative of the QoI wrt dof types
             */
            virtual void compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
            {
                MORIS_ERROR( false, "IQI::compute_dIQIdDof - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest wrt to dv types
             * @param[ in ] adQIdDv matrix to fill with derivative of the QI wrt dv types
             */
            virtual void compute_dQIdDv( Matrix< DDRMat > & adQIdDv )
            {
                MORIS_ERROR( false, "IQI::compute_dIQIdDv - This function does nothing. " );
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IQI_HPP_ */
