/*
 * cl_FEM_Stabilization_Parameter.hpp
 *
 *  Created on: Oct 18, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_STABILIZATION_PARAMETER_HPP_
#define SRC_FEM_CL_FEM_STABILIZATION_PARAMETER_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_Model.hpp"    //FEM/INT/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src

namespace moris
{
    namespace fem
    {
    class Cluster;

//------------------------------------------------------------------------------
        /**
         * Stabilization_Parameter
         */
        class Stabilization_Parameter
        {

        protected :

            // stabilization type
            fem::Penalty_Type mStabilizationType;

            // cluster pointer
            fem::Cluster * mCluster;

            // list of cluster measure enums
            moris::Cell< fem::Cluster_Measure > mClusterMeasures;

            // cluster measures
            // FIXME add enum for child class to select the needed ones
            real mMasterVolume;      // volume on master
            real mSlaveVolume;       // volume on slave
            real mInterfaceSurface;  // surface on master/slave interface
            real mElementSize;       // element size

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
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveDvTypes;

            // master and slave global dv type list
            moris::Cell< moris::Cell< MSI::Dv_Type > > mMasterGlobalDvTypes;
            moris::Cell< moris::Cell< MSI::Dv_Type > > mSlaveGlobalDvTypes;

            // master and slave global dv type maps
            Matrix< DDSMat > mMasterGlobalDvTypeMap;
            Matrix< DDSMat > mSlaveGlobalDvTypeMap;

            // master and slave dv field interpolators
            moris::Cell< Field_Interpolator* > mMasterDvFI;
            moris::Cell< Field_Interpolator* > mSlaveDvFI;

//            // master and slave property type lists
//            moris::Cell< fem::Property_Type > mMasterPropTypes;
//            moris::Cell< fem::Property_Type > mSlavePropTypes;

//            // master and slave global property type lists
//            moris::Cell< fem::Property_Type > mMasterGlobalPropTypes;
//            moris::Cell< fem::Property_Type > mSlaveGlobalPropTypes;

            // master and slave properties
            moris::Cell< std::shared_ptr< Property > > mMasterProp;
            moris::Cell< std::shared_ptr< Property > > mSlaveProp;

//            // master and slave constitutive type lists
//            moris::Cell< fem::Constitutive_Type > mMasterConstitutiveTypes;
//            moris::Cell< fem::Constitutive_Type > mSlaveConstitutiveTypes;

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
            /**
             * set the stabilization type
             * @param[ in ] aStabilizationType a penalty type
             */
            void set_penalty_type( fem::Penalty_Type aStabilizationType )
            {
                // set penalty type
                mStabilizationType = aStabilizationType;
            }

//------------------------------------------------------------------------------
            /**
             * get the stabilization type
             * @param[ out ] mStabilizationType a penalty type
             */
            fem::Penalty_Type get_stabilization_type()
            {
                // return stabilization type
                return mStabilizationType;
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
             */
            void reset_cluster_measures()
            {
                // FIXME cluster measures to reset volume, surface, ...
                mMasterVolume     = 0.5; // mCluster->compute_volume()
                mSlaveVolume      = 0.5; // mCluster->compute_volume()
                mInterfaceSurface = 1.0; // mCluster->compute_surface()
                mElementSize      = 1.0; // mCluster->compute_element_size( arg on how to compute )
            }

//------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags()
            {
                // reset the value flag
                mPPEval = true;

                // reset the master dof derivative flags
                uint tNumMasterDofTypes = mMasterGlobalDofTypes.size();
                mdPPdMasterDofEval.assign( tNumMasterDofTypes, true );

                // reset the slave dof derivative flags
                uint tNumSlaveDofTypes = mSlaveGlobalDofTypes.size();
                mdPPdSlaveDofEval.assign( tNumSlaveDofTypes, true );

                // reset the master dv derivative flags
                uint tNumMasterDvTypes = mMasterGlobalDvTypes.size();
                mdPPdMasterDvEval.assign( tNumMasterDvTypes, true );

                // reset the slave dv derivative flags
                uint tNumSlaveDvTypes = mSlaveGlobalDvTypes.size();
                mdPPdSlaveDvEval.assign( tNumSlaveDvTypes, true );

                // reset underlying master constitutive models
                for( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    tCM->reset_eval_flags();
                }

                // reset underlying slave constitutive models
                for( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    tCM->reset_eval_flags();
                }

                // reset underlying master properties
                for( std::shared_ptr< Property > tProp : mMasterProp )
                {
                    tProp->reset_eval_flags();
                }

                // reset underlying slave properties
                for( std::shared_ptr< Property > tProp : mSlaveProp )
                {
                    tProp->reset_eval_flags();
                }
            }

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
            void set_dv_type_list( const moris::Cell< moris::Cell< MSI::Dv_Type > > & aDvTypes,
                                    mtk::Master_Slave                                 aIsMaster = mtk::Master_Slave::MASTER )
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
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
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
                        MORIS_ASSERT( false, "Penalty_Parameter::get_global_dof_type_list - can only be master or slave." );
                        return mMasterGlobalDofTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * create a global dof type list including constitutive and property dependencies
             */
            void build_global_dof_type_list()
            {
                // MASTER-------------------------------------------------------
                // get the size of the dof type list
                uint tCounterMax = 0;

                // get number of dof types from penalty parameter
                tCounterMax += mMasterDofTypes.size();

                // get number of dof types from properties
                for ( std::shared_ptr< Property > tProperty : mMasterProp )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }

                // get number of dof types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    tCounterMax += tCM->get_global_dof_type_list().size();
                }

                // set size for the global dof type list
                mMasterGlobalDofTypes.resize( tCounterMax );

                // set a size for the checkList (used to avoid repeating a dof type)
                moris::Cell< sint > tCheckList( tCounterMax, -1 );

                // init total dof counter
                uint tCounter = 0;

                // get dof type from penalty parameter
                for ( uint iDOF = 0; iDOF < mMasterDofTypes.size(); iDOF++ )
                {
                    // put the dof type in the checklist
                    tCheckList( tCounter ) = static_cast< uint >( mMasterDofTypes( iDOF )( 0 ) );

                    // put the dof type in the global type list
                    mMasterGlobalDofTypes( tCounter ) = mMasterDofTypes( iDOF );

                    // update the dof counter
                    tCounter++;
                }

                // get dof type from properties
                for ( std::shared_ptr< Property > tProperty : mMasterProp )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            // put the dof type in the checklist
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                            // update dof counter
                            tCounter++;
                        }
                    }
                }

                // get dof type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            // put the dof type in the checklist
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                            // put the dof type in the global type list
                            mMasterGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                            // update dof counter
                            tCounter++;
                        }
                    }
                }

                // get the number of unique dof type groups for the penalty parameter
                mMasterGlobalDofTypes.resize( tCounter );

                // SLAVE--------------------------------------------------------
                // get the size of the dof type list
                tCounterMax = 0;

                // get number of dof types from penalty parameter
                tCounterMax += mSlaveDofTypes.size();

                // get number of dof types from properties
                for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                {
                    tCounterMax += tProperty->get_dof_type_list().size();
                }

                // get number of dof types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    tCounterMax += tCM->get_global_dof_type_list().size();
                }

                // set size for the global dof type list
                mSlaveGlobalDofTypes.resize( tCounterMax );

                // set a size for the checkList (used to avoid repeating a dof type)
                tCheckList.resize( tCounterMax, -1 );

                // init total dof counter
                tCounter = 0;

                // get dof type from penalty parameter
                for ( uint iDOF = 0; iDOF < mSlaveDofTypes.size(); iDOF++ )
                {
                    // put the dof type in the checklist
                    tCheckList( tCounter ) = static_cast< uint >( mSlaveDofTypes( iDOF )( 0 ) );

                    // put the dof type in the global type list
                    mSlaveGlobalDofTypes( tCounter ) = mSlaveDofTypes( iDOF );

                    // update the dof counter
                    tCounter++;
                }

                // get dof type from properties
                for ( std::shared_ptr< Property > tProperty : mSlaveProp )
                {
                    // get dof types for property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            // put the dof type in the checklist
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                            // update dof counter
                            tCounter++;
                        }
                    }
                }

                // get dof type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDOF = 0; iDOF < tActiveDofType.size(); iDOF++ )
                    {
                        // check enum is not already in the list
                        bool tCheck = false;
                        for( uint i = 0; i < tCounter; i++ )
                        {
                            tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tActiveDofType( iDOF )( 0 ) ) );
                        }

                        // if dof enum not in the list
                        if ( !tCheck )
                        {
                            // put the dof type in the checklist
                            tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );

                            // put the dof type in the global type list
                            mSlaveGlobalDofTypes( tCounter ) = tActiveDofType( iDOF );

                            // update dof counter
                            tCounter++;
                        }
                    }
                }

                // get the number of unique dof type groups for the penalty parameter
                mSlaveGlobalDofTypes.resize( tCounter );

                // number of global master and slave dof types
                uint tNumMasterGlobalDofTypes = mMasterGlobalDofTypes.size();
                uint tNumSlaveGlobalDofTypes  = mSlaveGlobalDofTypes.size();

                // set flag for evaluation
                mdPPdMasterDofEval.assign( tNumMasterGlobalDofTypes, true );
                mdPPdSlaveDofEval.assign( tNumSlaveGlobalDofTypes, true );

                // set storage for evaluation
                mdPPdMasterDof.resize( tNumMasterGlobalDofTypes );
                mdPPdSlaveDof.resize( tNumSlaveGlobalDofTypes );
            };

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
                        MORIS_ASSERT( false, "Penalty_Parameter::get_global_dof_type_map - can only be master or slave." );
                        return mMasterGlobalDofTypeMap;
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of master dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             */
            bool check_master_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // set bool for dependency
                bool tDofDependency = false;

                // get dof type index
                uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

                // if aDofType is an active dv type for the constitutive model
                if( tDofIndex < mMasterGlobalDofTypeMap.numel() && mMasterGlobalDofTypeMap( tDofIndex ) != -1 )
                {
                    // bool is set to true
                    tDofDependency = true;
                }
                // return bool for dependency
                return tDofDependency;
            }

//------------------------------------------------------------------------------
            /**
             * check dependency on a given group of slave dof types
             * @param[ in ]  aDofType       a group of dof types
             * @param[ out ] tDofDependency a bool true if dependency on dof type
             *
             */
            bool check_slave_dof_dependency( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // set bool for dependency
                bool tDofDependency = false;

                // get dof type index
                uint tDofIndex = static_cast< uint >( aDofType( 0 ) );

                // if aDofType is an active dv type for the constitutive model
                if( tDofIndex < mSlaveGlobalDofTypeMap.numel() && mSlaveGlobalDofTypeMap( tDofIndex ) != -1 )
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

                // set field interpolators for constitutive models
                for( std::shared_ptr< Constitutive_Model > tCM : this->get_constitutive_models( aIsMaster ) )
                {
                    // get the list of dof types for the CM
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tCMDofTypes = tCM->get_global_dof_type_list();

                    // get the number of dof type for the CM
                    uint tNumDofTypes = tCMDofTypes.size();

                    // set the size of the field interpolators list for the CM
                    moris::Cell< Field_Interpolator* > tCMFIs( tNumDofTypes, nullptr );

                    // loop over the dof types
                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
                    {
                        // get the dof type index in set
                        uint tDofIndexInSP = this->get_global_dof_type_map( aIsMaster )( static_cast< uint >( tCMDofTypes( iDof )( 0 ) ) );

                        // fill the field interpolators list for the CM
                        tCMFIs( iDof ) = this->get_dof_field_interpolators( aIsMaster )( tDofIndexInSP );
                    }

                    // set the field interpolators for the CM
                    tCM->set_dof_field_interpolators( tCMFIs );
                }

                // set field interpolators for properties
                for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
                {
                    // get the list of dof types for the property
                    moris::Cell< moris::Cell< MSI::Dof_Type > > tPropDofTypes = tProp->get_dof_type_list();

                    // get the number of dof type for the property
                    uint tNumDofTypes = tPropDofTypes.size();

                    // set the size of the field interpolators list for the property
                    moris::Cell< Field_Interpolator* > tPropFIs( tNumDofTypes, nullptr );

                    // loop over the dof types
                    for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
                    {
                        // get the dof type index in SP
                        uint tDofIndexInSP = this->get_global_dof_type_map( aIsMaster )( static_cast< uint >( tPropDofTypes( iDof )( 0 ) ) );

                        // fill the field interpolators list for the property
                        tPropFIs( iDof ) = this->get_dof_field_interpolators( aIsMaster )( tDofIndexInSP );
                    }

                    // set the field interpolators for the property
                    tProp->set_dof_field_interpolators( tPropFIs );
                }
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
                    tCM->set_geometry_interpolator( aGeometryInterpolator );
                }

                // set geometry interpolator for properties
                for( std::shared_ptr< Property > tProp : this->get_properties( aIsMaster ) )
                {
                    tProp->set_geometry_interpolator( aGeometryInterpolator );
                }
            }

//------------------------------------------------------------------------------
            void set_constitutive_models( moris::Cell< std::shared_ptr< Constitutive_Model > > aConstitutiveModels,
                                          mtk::Master_Slave                                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                this->get_constitutive_models( aIsMaster ) = aConstitutiveModels;
            }

//------------------------------------------------------------------------------
            moris::Cell< std::shared_ptr< Constitutive_Model > > & get_constitutive_models( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                return this->get_constitutive_models( aIsMaster );
            }

//------------------------------------------------------------------------------
            void set_properties( moris::Cell< std::shared_ptr< Property > > aProperties,
                                 mtk::Master_Slave                          aIsMaster = mtk::Master_Slave::MASTER )
            {
                this->get_properties( aIsMaster ) = aProperties;
            }

//------------------------------------------------------------------------------
            moris::Cell< std::shared_ptr< Property > > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                return this->get_properties( aIsMaster );
            }

//------------------------------------------------------------------------------
            /**
             * get global dv type list
             * @param[ out ] mGlobalDvTypes global list of dv type
             */
            const moris::Cell< moris::Cell< MSI::Dv_Type > > & get_global_dv_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
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
                        MORIS_ASSERT( false, "Penalty_Parameter::get_global_dv_type_list - can only be master or slave." );
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
                    tCounterMax += tProperty->get_dv_type_list().size();
                }

                // get number of dof types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    tCounterMax += tCM->get_global_dv_type_list().size();
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
                    // get dv types for property
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tProperty->get_dv_type_list();

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

                // get dof type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    // get dof types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tCM->get_global_dv_type_list();

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
                    tCounterMax += tProperty->get_dv_type_list().size();
                }

                // get number of dv types from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mSlaveCM )
                {
                    tCounterMax += tCM->get_global_dv_type_list().size();
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
                    // get dv types for property
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tProperty->get_dv_type_list();

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

                // get dv type from constitutive models
                for ( std::shared_ptr< Constitutive_Model > tCM : mMasterCM )
                {
                    // get dv types for constitutive model
                    moris::Cell< moris::Cell< MSI::Dv_Type > > tActiveDvType = tCM->get_global_dv_type_list();

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
            bool check_master_dv_dependency( const moris::Cell< MSI::Dv_Type > & aDvType )
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
            bool check_slave_dv_dependency( const moris::Cell< MSI::Dv_Type > & aDvType )
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
                    this->eval_PP();

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
            virtual void eval_PP()
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_PP - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter derivative wrt master dof
             * @param[ in ]  aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdPPdMasterDof penalty parameter derivative wrt master dof
             */
            const Matrix< DDRMat > & dPPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_master_dof_dependency( aDofType ), "Penalty_Parameter::dPPdMasterDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdPPdMasterDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dPPdMasterDOF( aDofType );

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
            virtual void eval_dPPdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_dPPdMasterDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter derivative wrt slave dof
             * @param[ in ]  aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdPPdSlaveDof penalty parameter derivative wrt master dof
             */
            const Matrix< DDRMat > & dPPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_slave_dof_dependency( aDofType ), "Penalty_Parameter::dPPdSlaveDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mSlaveGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdPPdSlaveDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dPPdSlaveDOF( aDofType );

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
            virtual void eval_dPPdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_dPPdSlaveDOF - This function does nothing. " );
            }

 //------------------------------------------------------------------------------
             /**
              * get the penalty parameter derivative wrt master dv
              * @param[ in ]  aDvTypes      a dv type wrt which the derivative is evaluated
              * @param[ out ] mdPPdMasterDv penalty parameter derivative wrt master dv
              */
             const Matrix< DDRMat > & dPPdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 // if aDofType is not an active dv type for the property
                 MORIS_ERROR( this->check_master_dv_dependency( aDvTypes ), "Penalty_Parameter::dPPdMasterDV - no dependency in this dv type." );

                 // get the dv index
                 uint tDvIndex = mMasterGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

                 // if the derivative has not been evaluated yet
                 if( mdPPdMasterDofEval( tDvIndex ) )
                 {
                     // evaluate the derivative
                     this->eval_dPPdMasterDV( aDvTypes );

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
             virtual void eval_dPPdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 MORIS_ERROR( false, " Penalty_Parameter::eval_dPPdMasterDV - This function does nothing. " );
             }

 //------------------------------------------------------------------------------
             /**
              * get the penalty parameter derivative wrt slave dv
              * @param[ in ]  aDvTypes     a dv type wrt which the derivative is evaluated
              * @param[ out ] mdPPdSlaveDv penalty parameter derivative wrt master dv
              */
             const Matrix< DDRMat > & dPPdSlaveDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 // if aDofType is not an active dv type for the property
                 MORIS_ERROR( this->check_slave_dv_dependency( aDvTypes ), "Penalty_Parameter::dPPdSlaveDV - no dependency in this dv type." );

                 // get the dv index
                 uint tDvIndex = mSlaveGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

                 // if the derivative has not been evaluated yet
                 if( mdPPdSlaveDvEval( tDvIndex ) )
                 {
                     // evaluate the derivative
                     this->eval_dPPdSlaveDV( aDvTypes );

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
             virtual void eval_dPPdSlaveDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 MORIS_ERROR( false, " Penalty_Parameter::eval_dPPdSlaveDV - This function does nothing. " );
             }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PENALTY_PARAMETER_HPP_ */
