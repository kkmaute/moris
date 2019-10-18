/*
 * cl_FEM_Penalty_Parameter.hpp
 *
 *  Created on: Oct 18, 2019
 *      Author: noel
 */
#ifndef SRC_FEM_CL_FEM_PENALTY_PARAMETER_HPP_
#define SRC_FEM_CL_FEM_PENALTY_PARAMETER_HPP_

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
//------------------------------------------------------------------------------
        /**
         * Penalty_Parameter
         */
        class Penalty_Parameter
        {

        protected :

            // cluster pointer
            fem::Cluster * mCluster;

            // cluster measures
            // FIXME add enum for child class to select the needed ones
            real mMasterVolume;      // volume on master
            real mSlaveVolume;       // volume on slave
            real mInterfaceSurface;  // surface on master/slave interface
            real mElementSize;       // element size

            // master and slave dof type lists
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // master and slave global dof type maps
            Matrix< DDSMat > mMasterGlobalDofTypeMap;
            Matrix< DDSMat > mSlaveGlobalDofTypeMap;

            // master and slave dof field interpolators
            moris::Cell< Field_Interpolator* > mMasterFI;
            moris::Cell< Field_Interpolator* > mSlaveFI;

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

            // master and slave property type lists
            moris::Cell< fem::Property_Type > mMasterPropTypes;
            moris::Cell< fem::Property_Type > mSlavePropTypes;

            // master and slave global property type lists
            moris::Cell< fem::Property_Type > mMasterGlobalPropTypes;
            moris::Cell< fem::Property_Type > mSlaveGlobalPropTypes;

            // master and slave properties
            moris::Cell< Property* > mMasterProp;
            moris::Cell< Property* > mSlaveProp;

            // master and slave constitutive type lists
            moris::Cell< fem::Constitutive_Type > mMasterConstitutiveTypes;
            moris::Cell< fem::Constitutive_Type > mSlaveConstitutiveTypes;

            // master and slave constitutive models
            moris::Cell< Constitutive_Model* > mMasterCM;
            moris::Cell< Constitutive_Model* > mSlaveCM;

            // flag for evaluation
            bool mPMEval = true;
            moris::Cell< bool > mdPMdMasterDofEval;
            moris::Cell< bool > mdPMdSlaveDofEval;
            moris::Cell< bool > mdPMdMasterDvEval;
            moris::Cell< bool > mdPMdSlaveDvEval;

            // storage
            Matrix< DDRMat > mPMVal = true;
            moris::Cell< Matrix< DDRMat > > mdPMdMasterDof;
            moris::Cell< Matrix< DDRMat > > mdPMdSlaveDof;
            moris::Cell< Matrix< DDRMat > > mdPMdMasterDv;
            moris::Cell< Matrix< DDRMat > > mdPMdSlaveDv;

//------------------------------------------------------------------------------
        public :

//------------------------------------------------------------------------------
            /**
             * constructor
             */
            Penalty_Parameter(){};

//------------------------------------------------------------------------------
            /**
             * virtual destructor
             */
            virtual ~Penalty_Parameter(){};

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
                mPMEval = true;

                // reset the master dof derivative flags
                uint tNumMasterDofTypes = mMasterGlobalDofTypes.size();
                mdPMdMasterDofEval.resize( tNumMasterDofTypes, true );

                // reset the slave dof derivative flags
                uint tNumSlaveDofTypes = mSlaveGlobalDofTypes.size();
                mdPMdSlaveDofEval.resize( tNumSlaveDofTypes, true );

                // reset the master dv derivative flags
                uint tNumMasterDvTypes = mMasterGlobalDvTypes.size();
                mdPMdMasterDvEval.resize( tNumMasterDvTypes, true );

                // reset the slave dv derivative flags
                uint tNumSlaveDvTypes = mSlaveGlobalDvTypes.size();
                mdPMdSlaveDvEval.resize( tNumSlaveDvTypes, true );
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
             * set property types
             * @param[ in ] aPropertyTypes a cell of property types
             * @param[ in ] aIsMaster enum for master or slave
             *
             */
            void set_property_type_list( const moris::Cell< fem::Property_Type  > & aPropertyTypes,
                                         mtk::Master_Slave                          aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterPropTypes = aPropertyTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlavePropTypes = aPropertyTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Penalty_Parameter::set_property_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * return a cell of property type
             * @param[ in ]  aIsMaster  enum master or slave
             * @param[ out ] mPropTypes a cell of property types
             */
            const moris::Cell< fem::Property_Type > & get_property_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master global property type list
                        return mMasterPropTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave global property type list
                        return mSlavePropTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Penalty_Parameter::get_property_type_list - can only be master or slave." );
                        return mMasterPropTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set constitutive types
             * @param[ in ] aConstitutiveTypes a cell of constitutive types
             * @param[ in ] aIsMaster enum for master or slave
             */
            void set_constitutive_type_list( const moris::Cell< fem::Constitutive_Type  > & aConstitutiveTypes,
                                             mtk::Master_Slave                              aIsMaster = mtk::Master_Slave::MASTER )
            {
                switch ( aIsMaster )
                {
                    case( mtk::Master_Slave::MASTER ) :
                    {
                        mMasterConstitutiveTypes = aConstitutiveTypes;
                        break;
                    }
                    case( mtk::Master_Slave::SLAVE ) :
                    {
                        mSlaveConstitutiveTypes = aConstitutiveTypes;
                        break;
                    }
                    default :
                    {
                        MORIS_ERROR( false, "Penalty_Parameter::set_constitutive_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * return a cell of constitutive types
             * @param[ in ]  aIsMaster enum master or slave
             * @param[ out ] mConstitutiveTypes a cell of constitutive types
             */
            const moris::Cell< fem::Constitutive_Type > & get_constitutive_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master constitutive type list
                        return mMasterConstitutiveTypes;
                        break;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave constitutive type list
                        return mSlaveConstitutiveTypes;
                        break;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "Penalty_Parameter::get_constitutive_type_list - can only be master or slave." );
                        return mMasterConstitutiveTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter value
             * @param[ out ] mPMVal penalty parameter value
             */
            const Matrix< DDRMat > & val()
            {
                // if the penalty parameter was not evaluated
                if( mPMEval )
                {
                    // evaluate the penalty parameter
                    this->eval_PM();
                }
                // return the penalty parameter value
                return mPMVal;
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter value
             */
            virtual void eval_PM()
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_PM - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter derivative wrt master dof
             * @param[ in ]  aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdPMdMasterDof penalty parameter derivative wrt master dof
             */
            const Matrix< DDRMat > & dPMdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_dof_dependency( aDofType ), "Penalty_Parameter::dPMdMasterDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mMasterGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdPMdMasterDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dPMdMasterDOF( aDofType );
                }

                // return the derivative
                return mdPMdMasterDof( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt master dof
             */
            virtual void eval_dPMdMasterDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_dPMdMasterDOF - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * get the penalty parameter derivative wrt slave dof
             * @param[ in ]  aDofTypes      a dof type wrt which the derivative is evaluated
             * @param[ out ] mdPMdSlaveDof penalty parameter derivative wrt master dof
             */
            const Matrix< DDRMat > & dPMdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofType )
            {
                // if aDofType is not an active dof type for the property
                MORIS_ERROR( this->check_dof_dependency( aDofType ), "Penalty_Parameter::dPMdSlaveDOF - no dependency in this dof type." );

                // get the dof index
                uint tDofIndex = mSlaveGlobalDofTypeMap( static_cast< uint >( aDofType( 0 ) ) );

                // if the derivative has not been evaluated yet
                if( mdPMdSlaveDofEval( tDofIndex ) )
                {
                    // evaluate the derivative
                    this->eval_dPMdSlaveDOF( aDofType );
                }

                // return the derivative
                return mdPMdSlaveDof( tDofIndex );
            }

//------------------------------------------------------------------------------
            /**
             * evaluate the penalty parameter derivative wrt slave dof
             */
            virtual void eval_dPMdSlaveDOF( const moris::Cell< MSI::Dof_Type > & aDofTypes )
            {
                MORIS_ERROR( false, " Penalty_Parameter::eval_dPMdSlaveDOF - This function does nothing. " );
            }

 //------------------------------------------------------------------------------
             /**
              * get the penalty parameter derivative wrt master dv
              * @param[ in ]  aDvTypes      a dv type wrt which the derivative is evaluated
              * @param[ out ] mdPMdMasterDv penalty parameter derivative wrt master dv
              */
             const Matrix< DDRMat > & dPMdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 // if aDofType is not an active dv type for the property
                 MORIS_ERROR( this->check_dof_dependency( aDvTypes ), "Penalty_Parameter::dPMdMasterDV - no dependency in this dv type." );

                 // get the dv index
                 uint tDvIndex = mMasterGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

                 // if the derivative has not been evaluated yet
                 if( mdPMdMasterDofEval( tDvIndex ) )
                 {
                     // evaluate the derivative
                     this->eval_dPMdMasterDV( aDvTypes );
                 }

                 // return the derivative
                 return mdPMdMasterDof( tDvIndex );
             }

 //------------------------------------------------------------------------------
             /**
              * evaluate the penalty parameter derivative wrt master dv
              */
             virtual void eval_dPMdMasterDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 MORIS_ERROR( false, " Penalty_Parameter::eval_dPMdMasterDV - This function does nothing. " );
             }

 //------------------------------------------------------------------------------
             /**
              * get the penalty parameter derivative wrt slave dv
              * @param[ in ]  aDvTypes     a dv type wrt which the derivative is evaluated
              * @param[ out ] mdPMdSlaveDv penalty parameter derivative wrt master dv
              */
             const Matrix< DDRMat > & dPMdSlaveDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 // if aDofType is not an active dof type for the property
                 MORIS_ERROR( this->check_dof_dependency( aDvTypes ), "Penalty_Parameter::dPMdSlaveDV - no dependency in this dv type." );

                 // get the dv index
                 uint tDvIndex = mSlaveGlobalDvTypeMap( static_cast< uint >( aDvTypes( 0 ) ) );

                 // if the derivative has not been evaluated yet
                 if( mdPMdSlaveDvEval( tDvIndex ) )
                 {
                     // evaluate the derivative
                     this->eval_dPMdSlaveDV( aDvTypes );
                 }

                 // return the derivative
                 return mdPMdSlaveDv( tDvIndex );
             }

 //------------------------------------------------------------------------------
             /**
              * evaluate the penalty parameter derivative wrt slave dv
              */
             virtual void eval_dPMdSlaveDV( const moris::Cell< MSI::Dv_Type > & aDvTypes )
             {
                 MORIS_ERROR( false, " Penalty_Parameter::eval_dPMdSlaveDV - This function does nothing. " );
             }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_PENALTY_PARAMETER_HPP_ */
