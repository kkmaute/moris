/*
 * cl_FEM_Stabilization_Parameter.hpp
 *
 *  Created on: Oct 18, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_STABILIZATION_PARAMETER_HPP_
#define SRC_FEM_CL_FEM_STABILIZATION_PARAMETER_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
//#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src

#include "cl_GEN_Dv_Enums.hpp"

namespace moris
{
    namespace fem
    {
    class Cluster;
    class Set;
    class Field_Interpolator_Manager;

//------------------------------------------------------------------------------
        /**
         * Stabilization_Parameter
         */
        class Stabilization_Parameter
        {
//------------------------------------------------------------------------------
        protected :

//------------------------------------------------------------------------------

            // fem set pointer
            fem::Set * mSet = nullptr;

            // field interpolator manager pointer
            Field_Interpolator_Manager * mMasterFIManager = nullptr;
            Field_Interpolator_Manager * mSlaveFIManager  = nullptr;

            // cluster pointer
            fem::Cluster * mCluster = nullptr;

            // list of parameters
            moris::Cell< Matrix< DDRMat > > mParameters;

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // bool for global dof type list and map
            bool mGlobalDofBuild = true;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave global dof type maps
            Matrix< DDSMat > mMasterGlobalDofTypeMap;
            Matrix< DDSMat > mSlaveGlobalDofTypeMap;

            // master and slave dof field interpolators
            moris::Cell< Field_Interpolator* > mMasterDofFI;
            moris::Cell< Field_Interpolator* > mSlaveDofFI;

            // master and slave dv type lists
            moris::Cell< moris::Cell< GEN_DV > > mMasterDvTypes;
            moris::Cell< moris::Cell< GEN_DV > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< GEN_DV > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< GEN_DV > > mSlaveGlobalDvTypes;

            // master and slave global dv type maps
            Matrix< DDSMat > mMasterGlobalDvTypeMap;
            Matrix< DDSMat > mSlaveGlobalDvTypeMap;

            // master and slave dv field interpolators
            moris::Cell< Field_Interpolator* > mMasterDvFI;
            moris::Cell< Field_Interpolator* > mSlaveDvFI;

            // master and slave properties
            moris::Cell< std::shared_ptr< Property > > mMasterProp;
            moris::Cell< std::shared_ptr< Property > > mSlaveProp;

            // master and slave constitutive models
            moris::Cell< std::shared_ptr< Constitutive_Model > > mMasterCM;
            moris::Cell< std::shared_ptr< Constitutive_Model > > mSlaveCM;

            // flag for evaluation
            bool mPPEval = true;
            moris::Cell< bool > mdPPdMasterDofEval;
            moris::Cell< bool > mdPPdSlaveDofEval;
            moris::Cell< bool > mdPPdMasterDvEval;
            moris::Cell< bool > mdPPdSlaveDvEval;

            // storage
            Matrix< DDRMat > mPPVal;
            moris::Cell< Matrix< DDRMat > > mdPPdMasterDof;
            moris::Cell< Matrix< DDRMat > > mdPPdSlaveDof;
            moris::Cell< Matrix< DDRMat > > mdPPdMasterDv;
            moris::Cell< Matrix< DDRMat > > mdPPdSlaveDv;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Stabilization_Parameter(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Stabilization_Parameter(){};

//------------------------------------------------------------------------------
            /*
             * set field interpolator manager pointer
             * @param[ in ] aFieldInteprolatorManager a field interpolator manager pointer
             */
            void set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager,
                                                 mtk::Master_Slave            aIsMaster = mtk::Master_Slave::MASTER );

//------------------------------------------------------------------------------
            /*
             * get field interpolator manager pointer
             * @param[ out ] aFieldInteprolatorManager a field interpolator manager pointer
             */
            Field_Interpolator_Manager * get_field_interpolator_manager( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case ( mtk::Master_Slave::MASTER ) :
                    {
                        return mMasterFIManager;
                    }

                    case ( mtk::Master_Slave::SLAVE ) :
                    {
                        return mSlaveFIManager;
                    }

                    default :
                    {
                        MORIS_ERROR( false, "Stabilization_Parameter::get_field_interpolator_manager - can only be master or slave");
                        return mMasterFIManager;
                    }
                }
            }

//------------------------------------------------------------------------------
            /*
             * set member set pointer
             * @param[ in ] aSetPointer a fem set pointer
             */
            void set_set_pointer( Set * aSetPointer )
            {
                mSet = aSetPointer;
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
            /**
             * set cluster
             * @param[ in ] aCluster a fem cluster pointer
             */
            void set_cluster( fem::Cluster * aCluster )
            {
                // set a cluster
                mCluster = aCluster;

                // reset cluster measures
                this->reset_cluster_measures();
            }

//------------------------------------------------------------------------------
            /**
             * reset cluster measures
             * NOTE: only implement if your stabilization parameter requires
             * cluster measure access. Otherwise no-op.
             */
            virtual
            void reset_cluster_measures()
            {

            }

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags();

//------------------------------------------------------------------------------
            /**
             * set dof types
             * @param[ in ] aDofTypes a list of group of dof types
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
                        MORIS_ERROR( false, "Penalty_Parameter::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dof types
             * @param[ in ] aIsMaster enum master or slave
             * @param[ out ] aDofTypes a list of group of dof types
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterDofTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveDofTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Penalty_Parameter::get_dof_type_list - can only be master or slave." );
                        return mMasterDofTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set dv types
             * @param[ in ] aDvTypes a list of group of dv types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_dv_type_list( const moris::Cell< moris::Cell< GEN_DV > > & aDvTypes,
                                    mtk::Master_Slave                           aIsMaster = mtk::Master_Slave::MASTER )
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
                        MORIS_ERROR( false, "Penalty_Parameter::set_dv_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of dv types
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] aDvTypes a list of group of dv types
             */
            const moris::Cell< moris::Cell< GEN_DV > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type list
                        return mMasterDvTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveDvTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Penalty_Parameter::get_dv_type_list - can only be master or slave." );
                        return mMasterDvTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * get global dof type list
             * @param[ out ] mGlobalDofTypes global list of dof type
             */
            const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                if( mGlobalDofBuild )
                {
                    // build the stabilization parameter global dof type list
                    this->build_global_dof_type_list();

                    // build the stabilization parameter global dof type map
                    this->build_global_dof_type_map();

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
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type list
                        return mSlaveGlobalDofTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dof_type_list - can only be master or slave." );
                        return mMasterGlobalDofTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * get a non unique list of dof type including
             * property, constitutive and stabilization dependencies
             * for both master and slave
             */
            void get_non_unique_dof_types( moris::Cell< MSI::Dof_Type >        & aDofTypes );
            void get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                  moris::Cell< GEN_DV >        & aDvTypes );
//------------------------------------------------------------------------------
            /**
             * create a global dof type list including constitutive and property dependencies
             */
            void build_global_dof_type_list();

//------------------------------------------------------------------------------
            /**
             * build global dof type map
             */
            void build_global_dof_type_map()
            {
                // MASTER-------------------------------------------------------
                // get number of global dof types
                uint tNumDofTypes = mMasterGlobalDofTypes.size();

                // determine the max Dof_Type enum
                sint tMaxEnum = 0;
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterGlobalDofTypes( iDOF )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dof_Type map size
                mMasterGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the Dof_Type map
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    // fill the property map
                    mMasterGlobalDofTypeMap( static_cast< int >( mMasterGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
                }

                // SLAVE-------------------------------------------------------
                // get number of global dof types
                tNumDofTypes = mSlaveGlobalDofTypes.size();

                // determine the max Dof_Type enum
                tMaxEnum = 0;
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveGlobalDofTypes( iDOF )( 0 ) ) );
                }
                tMaxEnum++;

                // set the dof type map size
                mSlaveGlobalDofTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the dof type map
                for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
                {
                    // fill the property map
                    mSlaveGlobalDofTypeMap( static_cast< int >( mSlaveGlobalDofTypes( iDOF )( 0 ) ), 0 ) = iDOF;
                }
            }

//------------------------------------------------------------------------------
            /**
             * get global dof type map
             */
            const Matrix< DDSMat > & get_global_dof_type_map( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dof type map
                        return mMasterGlobalDofTypeMap;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dof type map
                        return mSlaveGlobalDofTypeMap;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dof_type_map - can only be master or slave." );
                        return mMasterGlobalDofTypeMap;
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ in ]  aIsMaster      enum master or slave
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             */
             bool check_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType,
                                              mtk::Master_Slave              aIsMaster = mtk::Master_Slave::MASTER)
             {
                 // set bool for dependency
                 bool tDofDependency = false;

                 // get dof type index
                 uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

                 // if aDofType is an active dof type for the stabilization parameter
                 if( tDofIndex < this->get_global_dof_type_map( aIsMaster ).numel()
                     && this->get_global_dof_type_map( aIsMaster )( tDofIndex ) != -1 )
                 {
                     // bool is set to true
                     tDofDependency = true;
                 }
                 // return bool for dependency
                 return tDofDependency;
             }

//------------------------------------------------------------------------------
            /**
             * set field interpolators
             * @param[ in ] afieldInterpolators a lisy of field interpolator pointers
             * @param[ in ] aIsMaster           an enum for master or slave
             */
            void set_dof_field_interpolators( moris::Cell< Field_Interpolator* > aFieldInterpolators,
                                              mtk::Master_Slave                  aIsMaster = mtk::Master_Slave::MASTER )
            {
                // get input size
                uint tInputNumFI = aFieldInterpolators.size();

                // check input size
                MORIS_ASSERT( tInputNumFI == this->get_global_dof_type_list( aIsMaster ).size(),
                              "Stabilization_Parameter::set_dof_field_interpolators - wrong input size. " );

                // check dof field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < tInputNumFI; iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == this->get_global_dof_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "Stabilization_Parameter::set_dof_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_dof_field_interpolators( aIsMaster ) = aFieldInterpolators;

//                // set field interpolators for constitutive models
//                for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
//                {
//                    if( tCM != nullptr )
//                    {
//                        // get the list of dof types for the CM
//                        moris::Cell< moris::Cell< MSI::Dof_Type > > tCMDofTypes = tCM->get_global_dof_type_list();
//
//                        // get the number of dof type for the CM
//                        uint tNumDofTypes = tCMDofTypes.size();
//
//                        // set the size of the field interpolators list for the CM
//                        moris::Cell< Field_Interpolator* > tCMFIs( tNumDofTypes, nullptr );
//
//                        // loop over the dof types
//                        for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
//                        {
//                            // get the dof type index in set
//                            uint tDofIndexInSP = this->get_global_dof_type_map( aIsMaster )( static_cast< uint >( tCMDofTypes( iDof )( 0 ) ) );
//
//                            // fill the field interpolators list for the CM
//                            tCMFIs( iDof ) = this->get_dof_field_interpolators( aIsMaster )( tDofIndexInSP );
//                        }
//
//                        // set the field interpolators for the CM
//                        tCM->set_dof_field_interpolators( tCMFIs );
//                    }
//                }
//
//                // set field interpolators for properties
//                for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
//                {
//                    if( tProp != nullptr )
//                    {
//                        // get the list of dof types for the property
//                        moris::Cell< moris::Cell< MSI::Dof_Type > > tPropDofTypes = tProp->get_dof_type_list();
//
//                        // get the number of dof type for the property
//                        uint tNumDofTypes = tPropDofTypes.size();
//
//                        // set the size of the field interpolators list for the property
//                        moris::Cell< Field_Interpolator* > tPropFIs( tNumDofTypes, nullptr );
//
//                        // loop over the dof types
//                        for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
//                        {
//                            // get the dof type index in SP
//                            uint tDofIndexInSP = this->get_global_dof_type_map( aIsMaster )( static_cast< uint >( tPropDofTypes( iDof )( 0 ) ) );
//
//                            // fill the field interpolators list for the property
//                            tPropFIs( iDof ) = this->get_dof_field_interpolators( aIsMaster )( tDofIndexInSP );
//                        }
//
//                        // set the field interpolators for the property
//                        tProp->set_dof_field_interpolators( tPropFIs );
//                    }
//                }
            }

//------------------------------------------------------------------------------
            /**
             * get dof field interpolators
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aFieldInterpolators cell of dof field interpolator pointers
             */
            moris::Cell< Field_Interpolator* > & get_dof_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master field interpolator pointers
                        return mMasterDofFI;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave field interpolator pointers
                        return mSlaveDofFI;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Stabilization_Parameter::set_dof_field_interpolators - can only be master or slave." );
                        return mMasterDofFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * set geometry interpolator
             * @param[ in ] aGeometryInterpolator geometry interpolator pointers
             * @param[ in ] aIsMaster             enum for master or slave
             */
            void set_geometry_interpolator( Geometry_Interpolator* aGeometryInterpolator,
                                            mtk::Master_Slave      aIsMaster = mtk::Master_Slave::MASTER )
            {
                // set geometry interpolator for constitutive models
                for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
                {
                    if( tCM != nullptr )
                    {
                        tCM->set_geometry_interpolator( aGeometryInterpolator );
                    }
                }

                // set geometry interpolator for properties
                for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
                {
                    if( tProp != nullptr )
                    {
                        tProp->set_geometry_interpolator( aGeometryInterpolator );
                    }
                }
            }

//------------------------------------------------------------------------------
            virtual void set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                                                 std::string                           aConstitutiveString,
                                                 mtk::Master_Slave                     aIsMaster = mtk::Master_Slave::MASTER )
            {
                MORIS_ERROR( false, "Stabilization_Parameter::set_constitutive_model - This function does nothing." );
            }

//------------------------------------------------------------------------------
            moris::Cell< std::shared_ptr< Constitutive_Model > > & get_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master CM pointers
                        return mMasterCM;
                    }

                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave CM pointers
                        return mSlaveCM;
                    }

                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IWG::get_constitutive_models - can only be master or slave." );
                        return mMasterCM;
                    }
                }
            }

//------------------------------------------------------------------------------
            virtual void set_property( std::shared_ptr< Property > aProperty,
                                       std::string                 aPropertyString,
                                       mtk::Master_Slave           aIsMaster = mtk::Master_Slave::MASTER )
            {
                MORIS_ERROR( false, "Stabilization_Parameter::set_property - This function does nothing." );
            }

//------------------------------------------------------------------------------
            moris::Cell< std::shared_ptr< Property > > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
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
                        MORIS_ASSERT( false, "IWG::get_properties - can only be master or slave." );
                        return mMasterProp;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * get global dv type list
             * @param[ out ] mGlobalDvTypes global list of dv type
             */
            const moris::Cell< moris::Cell< GEN_DV > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global dv type list
                        return mMasterGlobalDvTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global dv type list
                        return mSlaveGlobalDvTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Stabilization_Parameter::get_global_dv_type_list - can only be master or slave." );
                        return mMasterGlobalDvTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * create a global dv type list including constitutive and property dependencies
             */
            void build_global_dv_type_list()
            {
                // MASTER-------------------------------------------------------
                // get the size of the dv type list
                uint tCounterMax = 0;

                // get number of dv types from penalty parameter
                tCounterMax += mMasterDvTypes.size();

                // get number of dv types from properties
                for ( std::shared_ptr< Property > tProperty : mMasterProp )
                {
                    if( tProperty != nullptr )
                    {
                        tCounterMax += tProperty->get_dv_type_list().size();
                    }
                }

                // get number of dof types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    if( tCM != nullptr )
                    {
                        tCounterMax += tCM->get_global_dv_type_list().size();
                    }
                }

                // set size for the global dv type list
                mMasterGlobalDvTypes.resize( tCounterMax );

                // set a size for the checkList (used to avoid repeating a dv type)
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dv counter
                uint tCounter = 0;

                // get dv type from penalty parameter
                for ( uint iDv = 0; iDv < mMasterDvTypes.size(); iDv++ )
                {
                    // put the dv type in the checklist
                    tCheckList( tCounter ) = static_cast< uint >( mMasterDvTypes( iDv )( 0 ) );

                    // put the dv type in the global type list
                    mMasterGlobalDvTypes( tCounter ) = mMasterDvTypes( iDv );

                    // update the dv counter
                    tCounter++;
                }

                // get dv type from properties
                for ( std::shared_ptr< Property > tProperty : mMasterProp )
                {
                    if( tProperty != nullptr )
                    {
                        // get dv types for property
                        moris::Cell< moris::Cell< GEN_DV > > tActiveDvType = tProperty->get_dv_type_list();

                        // loop on property dv type
                        for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                        {
                            // check enum is not already in the list
                            bool tCheck = false;
                            for( uint i = 0; i < tCounter; i++ )
                            {
                                tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                            }

                            // if dof enum not in the list
                            if ( !tCheck )
                            {
                                // put the dv type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                // put the dv type in the global type list
                                mMasterGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                // update dof counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get dof type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    if( tCM != nullptr )
                    {
                        // get dof types for constitutive model
                        moris::Cell< moris::Cell< GEN_DV > > tActiveDvType = tCM->get_global_dv_type_list();

                        // loop on property dv type
                        for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                        {
                            // check enum is not already in the list
                            bool tCheck = false;
                            for( uint i = 0; i < tCounter; i++ )
                            {
                                tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                            }

                            // if dv enum not in the list
                            if ( !tCheck )
                            {
                                // put the dv type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                // put the dv type in the global type list
                                mMasterGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                // update dv counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get the number of unique dv type groups for the penalty parameter
                mMasterGlobalDvTypes.resize( tCounter );

                // SLAVE--------------------------------------------------------
                // get the size of the dv type list
                tCounterMax = 0;

                // get number of dv types from penalty parameter
                tCounterMax += mSlaveDvTypes.size();

                // get number of dv types from properties
                for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                {
                    if( tProperty != nullptr )
                    {
                        tCounterMax += tProperty->get_dv_type_list().size();
                    }
                }

                // get number of dv types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    if( tCM != nullptr )
                    {
                        tCounterMax += tCM->get_global_dv_type_list().size();
                    }
                }

                // set size for the global dv type list
                mSlaveGlobalDvTypes.resize( tCounterMax );

                // set a size for the checkList (used to avoid repeating a dv type)
                tCheckList.resize( tCounterMax, -1 );

                // init total dv counter
                tCounter = 0;

                // get dv type from penalty parameter
                for ( uint iDv = 0; iDv < mSlaveDvTypes.size(); iDv++ )
                {
                    // put the dv type in the checklist
                    tCheckList( tCounter ) = static_cast< uint >( mSlaveDvTypes( iDv )( 0 ) );

                    // put the dv type in the global type list
                    mSlaveGlobalDvTypes( tCounter ) = mSlaveDvTypes( iDv );

                    // update the dv counter
                    tCounter++;
                }

                // get dv type from properties
                for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                {
                    if( tProperty != nullptr )
                    {
                        // get dv types for property
                        moris::Cell< moris::Cell< GEN_DV > > tActiveDvType = tProperty->get_dv_type_list();

                        // loop on property dv type
                        for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                        {
                            // check enum is not already in the list
                            bool tCheck = false;
                            for( uint i = 0; i < tCounter; i++ )
                            {
                                tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                            }

                            // if dv enum not in the list
                            if ( !tCheck )
                            {
                                // put the dv type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                // put the dv type in the global type list
                                mSlaveGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                // update dv counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get dv type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    if( tCM != nullptr )
                    {
                        // get dv types for constitutive model
                        moris::Cell< moris::Cell< GEN_DV > > tActiveDvType = tCM->get_global_dv_type_list();

                        // loop on property dv type
                        for ( uint iDv = 0; iDv < tActiveDvType.size(); iDv++ )
                        {
                            // check enum is not already in the list
                            bool tCheck = false;
                            for( uint i = 0; i < tCounter; i++ )
                            {
                                tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDvType( iDv )( 0 ) ) );
                            }

                            // if dv enum not in the list
                            if ( !tCheck )
                            {
                                // put the dv type in the checklist
                                tCheckList( tCounter ) = static_cast< uint >( tActiveDvType( iDv )( 0 ) );

                                // put the dv type in the global type list
                                mSlaveGlobalDvTypes( tCounter ) = tActiveDvType( iDv );

                                // update dv counter
                                tCounter++;
                            }
                        }
                    }
                }

                // get the number of unique dv type groups for the penalty parameter
                mSlaveGlobalDvTypes.resize( tCounter );

                // build global dv type map
                this->build_global_dv_type_map();

                // number of global master and slave dv types
                uint tNumMasterGlobalDvTypes = mMasterGlobalDvTypes.size();
                uint tNumSlaveGlobalDvTypes  = mSlaveGlobalDvTypes.size();

                // set flag for evaluation
                mdPPdMasterDvEval.assign( tNumMasterGlobalDvTypes, true );
                mdPPdSlaveDvEval.assign( tNumSlaveGlobalDvTypes, true );

                // set storage for evaluation
                mdPPdMasterDv.resize( tNumMasterGlobalDvTypes );
                mdPPdSlaveDv.resize( tNumSlaveGlobalDvTypes );
            };

//------------------------------------------------------------------------------
            /**
             * build global dv type map
             */
            void build_global_dv_type_map()
            {
                // MASTER-------------------------------------------------------
                // get number of global dof types
                uint tNumDvTypes = mMasterGlobalDvTypes.size();

                // determine the max Dv_Type enum
                sint tMaxEnum = 0;
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mMasterGlobalDvTypes( iDv )( 0 ) ) );
                }
                tMaxEnum++;

                // set the Dv_Type map size
                mMasterGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the Dv_Type map
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    // fill the property map
                    mMasterGlobalDvTypeMap( static_cast< int >( mMasterGlobalDvTypes( iDv )( 0 ) ), 0 ) = iDv;
                }

                // SLAVE-------------------------------------------------------
                // get number of global dv types
                tNumDvTypes = mSlaveGlobalDvTypes.size();

                // determine the max Dv_Type enum
                tMaxEnum = 0;
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    tMaxEnum = std::max( tMaxEnum, static_cast< int >( mSlaveGlobalDvTypes( iDv )( 0 ) ) );
                }
                tMaxEnum++;

                // set the dv type map size
                mSlaveGlobalDvTypeMap.set_size( tMaxEnum, 1, -1 );

                // fill the dv type map
                for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
                {
                    // fill the property map
                    mSlaveGlobalDvTypeMap( static_cast< int >( mSlaveGlobalDvTypes( iDv )( 0 ) ), 0 ) = iDv;
                }
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of master dv types
             * @param[ in ]  aDvType       a group of dv types
             * @param[ out ] tDvDependency a bool true if dependency on dv type
             */
            bool check_master_dv_dependency( const moris::Cell< GEN_DV > & aDvType )
            {
                // set bool for dependency
                bool tDvDependency = false;

                // get dv type index
                uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

                // if aDvType is an active dv type for the constitutive model
                if( tDvIndex < mMasterGlobalDvTypeMap.numel() && mMasterGlobalDvTypeMap( tDvIndex ) != -1 )
                {
                    // bool is set to true
                    tDvDependency = true;
                }
                // return bool for dependency
                return tDvDependency;
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of slave dv types
             * @param[ in ]  aDvType       a group of dv types
             * @param[ out ] tDvDependency a bool true if dependency on dv type
             *
             */
            bool check_slave_dv_dependency( const moris::Cell< GEN_DV > & aDvType )
            {
                // set bool for dependency
                bool tDvDependency = false;

                // get dv type index
                uint tDvIndex = static_cast< uint >( aDvType( 0 ) );

                // if aDvType is an active dv type for the constitutive model
                if( tDvIndex < mSlaveGlobalDvTypeMap.numel() && mSlaveGlobalDvTypeMap( tDvIndex ) != -1 )
                {
                    // bool is set to true
                    tDvDependency = true;
                }
                // return bool for dependency
                return tDvDependency;
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter value
             * @param[ out ] mPPVal penalty parameter value
             */
            const Matrix< DDRMat > & val()
            {
                // if the penalty parameter was not evaluated
                if( mPPEval )
                {
                    // evaluate the penalty parameter
                    this->eval_SP();

                    // set bool for evaluation
                    mPPEval = false;
                }
                // return the penalty parameter value
                return mPPVal;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter value
             */
            virtual void eval_SP()
            {
                MORIS_ERROR( false, " Stabilization_Parameter::eval_PP - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter derivative wrt master dof
             * @param[ in ]  aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdPPdMasterDof penalty parameter derivative wrt master dof
             */
            const Matrix< DDRMat > & dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_dof_dependency( aDofType, mtk::Master_Slave::MASTER ),
                             "Stabilization_Parameter::dPPdMasterDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdPPdMasterDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dSPdMasterDOF( aDofType );

                    // set bool for evaluation
                    mdPPdMasterDofEval( tDofIndex ) = false;
                }

                // return the derivative
                return mdPPdMasterDof( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt master dof
             */
            virtual void eval_dSPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_dSPdMasterDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter derivative wrt slave dof
             * @param[ in ]  aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdPPdSlaveDof penalty parameter derivative wrt master dof
             */
            const Matrix< DDRMat > & dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_dof_dependency( aDofType, mtk::Master_Slave::SLAVE ),
                             "Stabilization_Parameter::dSPdSlaveDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mSlaveGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdPPdSlaveDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dSPdSlaveDOF( aDofType );

                    // set bool for evaluation
                    mdPPdSlaveDofEval( tDofIndex ) = false;
                }

                // return the derivative
                return mdPPdSlaveDof( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt slave dof
             */
            virtual void eval_dSPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Stabilization_Parameter::eval_dSPdSlaveDOF - This function does nothing. " );
            }

 //------------------------------------------------------------------------------
             /**
              * get the penalty parameter derivative wrt master dv
              * @param[ in ]  aDvTypes      a dv type wrt which the derivative is evaluated
              * @param[ out ] mdPPdMasterDv penalty parameter derivative wrt master dv
              */
             const Matrix< DDRMat > & dSPdMasterDV( const moris::Cell< GEN_DV > & aDvTypes )
             {
                 // if aDofType is not an active dv type for the property
                 MORIS_ERROR( this->check_master_dv_dependency( aDvTypes ), "Penalty_Parameter::dPPdMasterDV - no dependency in this dv type." );

                 // get the dv index
                 uint tDvIndex = mMasterGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

                 // if the derivative has not been evaluated yet
                 if( mdPPdMasterDofEval( tDvIndex ) )
                 {
                     // evaluate the derivative
                     this->eval_dSPdMasterDV( aDvTypes );

                     // set bool for evaluation
                     mdPPdMasterDofEval( tDvIndex ) = false;
                 }

                 // return the derivative
                 return mdPPdMasterDof( tDvIndex );
             }

 //------------------------------------------------------------------------------
             /**
              * evaluate the penalty parameter derivative wrt master dv
              */
             virtual void eval_dSPdMasterDV( const moris::Cell< GEN_DV > & aDvTypes )
             {
                 MORIS_ERROR( false, " Stabilization_Parameter::eval_dSPdMasterDV - This function does nothing. " );
             }

 //------------------------------------------------------------------------------
             /**
              * get the penalty parameter derivative wrt slave dv
              * @param[ in ]  aDvTypes     a dv type wrt which the derivative is evaluated
              * @param[ out ] mdPPdSlaveDv penalty parameter derivative wrt master dv
              */
             const Matrix< DDRMat > & dSPdSlaveDV( const moris::Cell< GEN_DV > & aDvTypes )
             {
                 // if aDofType is not an active dv type for the property
                 MORIS_ERROR( this->check_slave_dv_dependency( aDvTypes ),
                              "Stabilization_Parameter::dSPdSlaveDV - no dependency in this dv type." );

                 // get the dv index
                 uint tDvIndex = mSlaveGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

                 // if the derivative has not been evaluated yet
                 if( mdPPdSlaveDvEval( tDvIndex ) )
                 {
                     // evaluate the derivative
                     this->eval_dSPdSlaveDV( aDvTypes );

                     // set bool for evaluation
                     mdPPdSlaveDvEval( tDvIndex ) = false;
                 }

                 // return the derivative
                 return mdPPdSlaveDv( tDvIndex );
             }

 //------------------------------------------------------------------------------
             /**
              * evaluate the penalty parameter derivative wrt slave dv
              */
             virtual void eval_dSPdSlaveDV( const moris::Cell< GEN_DV > & aDvTypes )
             {
                 MORIS_ERROR( false, " Stabilization_Parameter::eval_dSPdSlaveDV - This function does nothing. " );
             }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_STABILIZATION_PARAMETER_HPP_ */
