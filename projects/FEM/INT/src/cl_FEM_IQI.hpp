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
#include "cl_GEN_Dv_Enums.hpp"

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
            moris::Cell< moris::Cell< GEN_DV > > mMasterDvTypes;
            moris::Cell< moris::Cell< GEN_DV > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< GEN_DV > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< GEN_DV > > mSlaveGlobalDvTypes;

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

            std::string mName;

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
            void print_names()
            {
                std::cout<<"----------"<<std::endl;
                std::cout<<"IQI: "<<mName<<std::endl;

                // properties
                for( uint iProp = 0; iProp < mMasterProp.size(); iProp++ )
                {
                    if( mMasterProp( iProp ) != nullptr )
                    {
                        std::cout<<"Master property: "<<mMasterProp( iProp )->get_name()<<std::endl;
                    }
                }
                for( uint iProp = 0; iProp < mSlaveProp.size(); iProp++ )
                {
                    if( mSlaveProp( iProp ) != nullptr )
                    {
                        std::cout<<"Slave property:  "<<mSlaveProp( iProp )->get_name()<<std::endl;
                    }
                }

                // CM
                for( uint iCM = 0; iCM < mMasterCM.size(); iCM++ )
                {
                    if( mMasterCM( iCM ) != nullptr )
                    {
                        std::cout<<"Master CM: "<<mMasterCM( iCM )->get_name()<<std::endl;
                    }
                }
                for( uint iCM = 0; iCM < mSlaveCM.size(); iCM++ )
                {
                    if( mSlaveCM( iCM ) != nullptr )
                    {
                        std::cout<<"Slave CM:  "<<mSlaveCM( iCM )->get_name()<<std::endl;
                    }
                }

                // SP
                for( uint iSP = 0; iSP < mStabilizationParam.size(); iSP++ )
                {
                    if( mStabilizationParam( iSP ) != nullptr )
                    {
                        std::cout<<"SP: "<<mStabilizationParam( iSP )->get_name()<<std::endl;
                    }
                }
                std::cout<<"----------"<<std::endl;
            }

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
            void get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                  moris::Cell< GEN_DV >  & aDvTypes );

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
                    this->build_global_dof_and_dv_type_lists();

                    // update build flag
                    mGlobalDofBuild = false;
                    mGlobalDvBuild  = false;
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
            void set_dv_type_list( const moris::Cell< moris::Cell< GEN_DV > > & aDvTypes,
                                   mtk::Master_Slave                            aIsMaster = mtk::Master_Slave::MASTER )
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
            moris::Cell< moris::Cell< GEN_DV > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
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
            void get_non_unique_global_dv_type_list( moris::Cell< GEN_DV > & aGlobalDvTypeList );

//------------------------------------------------------------------------------
            /**
             * get a unique global dv type list including
             * IQI, property, constitutive and stabilization dependencies
             * @param[ in ] aIsMaster enum master or slave
             */
            moris::Cell< moris::Cell< GEN_DV > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // if the global list was not yet built
                if( mGlobalDvBuild )
                {
                    // build global dv type list
                    this->build_global_dof_and_dv_type_lists();

                    // update build flag
                    mGlobalDvBuild = false;
                    mGlobalDofBuild = false;
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
                MORIS_ERROR( false, "IQI::compute_QI - Not implemented for base class. " );
            }

//------------------------------------------------------------------------------
            /**
             * compute the quantities of interest
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_QI( real aWStar )
            {
                MORIS_ERROR( false, "IQI::compute_QI - Not implemented for base class. " );
            }

//------------------------------------------------------------------------------
            /**
             * get requested dof type
             * @param[ in ] mRequestedDofType list of requested dof type
             */
            moris::Cell < enum MSI::Dof_Type > get_requested_dof_types();

//------------------------------------------------------------------------------
            /**
             * get requested dof type
             * @param[ in ] mRequestedMasterDofTypes list of requested master dof type
             * @param[ in ] mRequestedSlaveDofTypes list of requested master dof type
             */
            void build_requested_dof_type_lists();

//------------------------------------------------------------------------------
            /**
             * compute the derivative of the quantities of interest
             * wrt requested dof types
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            virtual void compute_dQIdu( moris::real aWStar )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdu - Not implemented for base class. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest
             * wrt to requested dof types by finite difference
             * @param[ in ] adQIdDuFD     matrix to fill with derivatives of the QI wrt dof types
             * @param[ in ] aWStar        weight associated to the evaluation point
             * @param[ in ] aPerturbation real for relative perturbation of the dof values
             */
            void compute_dQIdu_FD( Matrix< DDRMat > & adQIdDuFD,
                                   real               aWStar,
                                   real               aPerturbation );

//------------------------------------------------------------------------------
            /**
             * check the derivative of the quantity of interest wrt to dof types
             * with evaluation by finite difference
             * @param[ in ] aPerturbation real for perturbation of the dof values
             * @param[ in ] aEpsilon      real for tolerance
             * @param[ in ] adQIdu        matrix to fill with derivative of QI wrt dof types
             * @param[ in ] adQIduFD      matrix to fill with derivative of QI wrt dof types
             *                            evaluated by finite difference
             */
            bool check_dQIdu_FD( real               aWStar,
                                 real               aPerturbation,
                                 real               aEpsilon,
                                 Matrix< DDRMat > & adQIdu,
                                 Matrix< DDRMat > & adQIduFD );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest wrt to dv types
             * @param[ in ] adQIdp matrix to fill with derivative of the QI wrt dv
             */
            virtual void compute_dQIdp( Matrix< DDRMat > & adQIdp )
            {
                MORIS_ERROR( false, "IQI::compute_dQIdp - Not implemented for base class. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest
             * wrt to material dv by finite difference
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation dv relative perturbation
             * @param[ in ] adQIdpMatFD   cell of matrix for dQIdpMat to fill
             */
            void compute_dQIdp_FD_material( moris::real        aWStar,
                                            moris::real        aPerturbation,
                                            Matrix< DDRMat > & adQIdpMatFD );

//------------------------------------------------------------------------------
            /**
             * evaluate the derivative of the quantity of interest
             * wrt to geometry dv by finite difference
             * @param[ in ] aWStar        weight associated to evaluation point
             * @param[ in ] aPerturbation dv relative perturbation
             * @param[ in ] aIsActive     cell of vectors for active dv
             * @param[ in ] adQIdpGeoFD   cell of matrix for dRdpGeo to fill
             */
            void compute_dQIdp_FD_geometry( moris::real                       aWStar,
                                            moris::real                       aPerturbation,
                                            moris::Cell< Matrix< DDSMat > > & aIsActive,
                                            Matrix< DDRMat >                & adQIdpGeoFD );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IQI_HPP_ */
