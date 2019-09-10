/*
 * cl_FEM_IWG.hpp
 *
 *  Created on: Aug 13, 2018
 *      Author: messe
 */
#ifndef SRC_FEM_CL_FEM_IWG_HPP_
#define SRC_FEM_CL_FEM_IWG_HPP_

#include "typedefs.hpp"                     //MRS/COR/src
#include "linalg_typedefs.hpp"              //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src
#include "cl_Matrix.hpp"                    //LNA/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Property.hpp"    //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"        //FEM/MSI/src
#include "cl_FEM_Enums.hpp"        //FEM/MSI/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        /**
         * Integrand of Weak Form of Governing Equations
         */
        class IWG
        {
        protected :

            // nodal weak BCs
            Matrix< DDRMat > mNodalWeakBCs;

            // normal
            Matrix< DDRMat > mNormal;

            // residual dof type
            moris::Cell< MSI::Dof_Type > mResidualDofType;

            // master and slave dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveDofTypes;

            // master and slave field interpolators
            moris::Cell< Field_Interpolator* > mMasterFI;
            moris::Cell< Field_Interpolator* > mSlaveFI;

            // master and slave property type list
            moris::Cell< fem::Property_Type > mMasterPropTypes;
            moris::Cell< fem::Property_Type > mSlavePropTypes;

            // master and slave properties
            moris::Cell< Property* > mMasterProp;
            moris::Cell< Property* > mSlaveProp;

            // master and slave global dof type list
            moris::Cell< moris::Cell< MSI::Dof_Type > > mMasterGlobalDofTypes;
            moris::Cell< moris::Cell< MSI::Dof_Type > > mSlaveGlobalDofTypes;

            // FIXME temporary until other way
            uint mSpaceDim = 3;

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
            void set_normal( Matrix< DDRMat > & aNormal )
            {
                mNormal = aNormal;
            }

//------------------------------------------------------------------------------
            /**
             * set residual dof type
             * @param[ in ] aResidualdofType a cell of residual dof types
             */
            void set_residual_dof_type( moris::Cell< MSI::Dof_Type > & aResidualDofType )
            {
                mResidualDofType = aResidualDofType;
            }

//------------------------------------------------------------------------------
            /**
             * returns a dof type for the residual
             */
            const moris::Cell< MSI::Dof_Type > & get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * set IWG active dof types
             * @param[ in ] aDofTypes a cell of cell of dof types
             * @param[ in ] aIsMaster enum for master or slave
             *
             */
            void set_dof_type_list( moris::Cell< moris::Cell< MSI::Dof_Type > > & aDofTypes,
                                    mtk::Master_Slave                             aIsMaster = mtk::Master_Slave::MASTER)
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
                        MORIS_ERROR( false, "IWG::set_dof_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * returns a cell of dof types active for the IWG
             * @param[ in ] aIsMaster enum master or slave
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
                        MORIS_ASSERT( false, "IWG::get_dof_type_list - can only be master or slave." );
                        return mMasterDofTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set IWG active property types
             * @param[ in ] aPropertyTypes a cell of cell of property types
             * @param[ in ] aIsMaster enum for master or slave
             *
             */
            void set_property_type_list( moris::Cell< fem::Property_Type  > & aPropertyTypes,
                                         mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER)
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
                        MORIS_ERROR( false, "IWG::set_property_type_list - can only be MASTER or SLAVE.");
                        break;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * returns a cell of property type active for the IWG
             * @param[ in ] aIsMaster enum master or slave
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
                        MORIS_ASSERT( false, "IWG::get_property_type_list - can only be master or slave." );
                        return mMasterPropTypes;
                        break;
                    }
                }
            };

//------------------------------------------------------------------------------
            /**
             * set field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             * @param[ in ] aIsMaster           enum for master or slave
             */
            void set_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators,
                                          mtk::Master_Slave                    aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == this->get_global_dof_type_list( aIsMaster ).size(),
                              "IWG::set_field_interpolators - wrong input size. " );

                // check field interpolator type
                bool tCheckFI = true;
                for( uint iFI = 0; iFI < aFieldInterpolators.size(); iFI++ )
                {
                    tCheckFI = tCheckFI && ( aFieldInterpolators( iFI )->get_dof_type()( 0 ) == this->get_global_dof_type_list( aIsMaster )( iFI )( 0 ) );
                }
                MORIS_ASSERT( tCheckFI, "IWG::set_field_interpolators - wrong field interpolator dof type. ");

                // set field interpolators
                this->get_field_interpolators( aIsMaster ) = aFieldInterpolators;
            }

//------------------------------------------------------------------------------
            /**
             * get field interpolators
             * @param[ in ]  aIsMaster           enum master or slave
             * @param[ out ] aFieldInterpolators cell of field interpolator pointers
             */
            moris::Cell< Field_Interpolator* > & get_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // switch on master/slave
                switch( aIsMaster )
                {
                    // if master
                    case( mtk::Master_Slave::MASTER ):
                    {
                        // return master field interpolator pointers
                        return mMasterFI;
                    }
                    // if slave
                    case( mtk::Master_Slave::SLAVE ):
                    {
                        // return slave field interpolator pointers
                        return mSlaveFI;
                    }
                    // if none
                    default:
                    {
                        MORIS_ASSERT( false, "IWG::set_field_interpolators - can only be master or slave." );
                        return mMasterFI;
                    }
                }
            }

//------------------------------------------------------------------------------
            /**
             * check that field interpolators were assigned
             * @param[ in ]  aIsMaster enum master or slave
             */
             void check_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check field interpolators cell size
                 MORIS_ASSERT( this->get_field_interpolators( aIsMaster ).size() == this->get_dof_type_list( aIsMaster ).size(),
                               "IWG::check_field_interpolators - master, wrong FI size. " );

                // loop over the field interpolator pointers
                for( uint iFI = 0; iFI < this->get_dof_type_list( aIsMaster ).size(); iFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT( this->get_field_interpolators( aIsMaster )( iFI ) != nullptr,
                                  "IWG::check_field_interpolators - master, FI missing. " );
                }
             }

//------------------------------------------------------------------------------
             /**
              * set properties
              * @param[ in ] aProperties cell of property pointers
              * @param[ in ] aIsMaster   enum master or slave
              */
             void set_properties( moris::Cell< Property* > & aProperties,
                                  mtk::Master_Slave          aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check input size
                 MORIS_ASSERT( aProperties.size() == this->get_property_type_list( aIsMaster ).size(),
                               "IWG::set_properties - master, wrong input size. " );

                 // check property type
                 bool tCheckProp = true;
                 for( uint iProp = 0; iProp < aProperties.size(); iProp++ )
                 {
                     tCheckProp = tCheckProp && ( aProperties( iProp )->get_property_type() == this->get_property_type_list( aIsMaster)( iProp ) );
                 }
                 MORIS_ASSERT( tCheckProp, "IWG::set_properties - wrong property type. ");

                 // set properties
                 this->get_properties( aIsMaster ) = aProperties;

                 // create a global dof type list
                 this->build_global_dof_type_list( aIsMaster );
             }

//------------------------------------------------------------------------------
             /**
              * get properties
              * @param[ in ]  aIsMaster   enum master or slave
              * @param[ out ] aProperties cell of property pointers
              */
             moris::Cell< Property* > & get_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // switch on master/slave
                 switch( aIsMaster )
                 {
                     // if master
                     case( mtk::Master_Slave::MASTER ):
                     {
                         // return master property pointers
                         return mMasterProp;
                         break;
                     }
                     // if slave
                     case( mtk::Master_Slave::SLAVE ):
                     {
                         // return slave property pointers
                         return mSlaveProp;
                         break;
                     }
                     // if none
                     default:
                     {
                         MORIS_ASSERT( false, "IWG::get_properties - can only be master or slave." );
                         return mMasterProp;
                         break;
                     }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * check that properties were assigned
              * @param[ in ]  aIsMaster enum master or slave
              */
              void check_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // check master properties cell size
                  MORIS_ASSERT( this->get_properties( aIsMaster ).size() == this->get_property_type_list( aIsMaster ).size(),
                                "IWG::check_properties - wrong property size. " );

                  // loop over all master properties and check that they are assigned
                  for( uint iProp = 0; iProp < mMasterPropTypes.size(); iProp++ )
                  {
                      MORIS_ASSERT( this->get_properties( aIsMaster )( iProp ) != nullptr,
                                    "IWG::check_properties - property missing. " );
                  }
              }

//------------------------------------------------------------------------------
              /**
               * create a global dof type list including IWG and properties dependencies
               * @param[ in ] aIsMaster enum master or slave
               */
              void build_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // set the size of the dof type list
                  uint tCounterMax = this->get_dof_type_list( aIsMaster ).size();

                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      tCounterMax += tProperty->get_dof_type_list().size();
                  }
                  this->get_global_dof_type_list( aIsMaster ).resize( tCounterMax );
                  moris::Cell< sint > tCheckList( tCounterMax, -1 );

                  // init total dof counter
                  uint tCounter = 0;

                  // get active dof type for IWG
                  for ( uint iDOF = 0; iDOF < this->get_dof_type_list( aIsMaster ).size(); iDOF++ )
                  {
                      tCheckList( tCounter ) = static_cast< uint >( this->get_dof_type_list( aIsMaster )( iDOF )( 0 ) );
                      this->get_global_dof_type_list( aIsMaster )( tCounter ) = this->get_dof_type_list( aIsMaster )( iDOF );
                      tCounter++;
                  }

                  for ( Property* tProperty : this->get_properties( aIsMaster ) )
                  {
                      // get active dof type
                      moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_dof_type_list();

                      for ( uint iDOF = 0; iDOF < tProperty->get_dof_type_list().size(); iDOF++ )
                      {
                          // check enum is not already in the list
                          bool tCheck = false;
                          for( uint i = 0; i < tCounter; i++ )
                          {
                              tCheck = tCheck || equal_to( tCheckList( i ), static_cast< uint >( tProperty->get_dof_type_list()( iDOF )( 0 ) ) );
                          }

                          // if dof enum not in the list
                          if ( !tCheck )
                          {
                              tCheckList( tCounter ) = static_cast< uint >( tProperty->get_dof_type_list()( iDOF )( 0 ) );
                              this->get_global_dof_type_list( aIsMaster )( tCounter ) = tProperty->get_dof_type_list()( iDOF );
                              tCounter++;
                          }
                      }
                  }

                  // get the number of unique dof type groups, i.e. the number of interpolators
                  this->get_global_dof_type_list( aIsMaster ).resize( tCounter );
              };

//------------------------------------------------------------------------------
              /**
               * get global dof type list
               * @param[ in ] aIsMaster enum master or slave
               */
              moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_type_list( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
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
                          MORIS_ASSERT( false, "IWG::get_global_dof_type_list - can only be master or slave." );
                          return mMasterGlobalDofTypes;
                          break;
                      }
                  }
              };

//------------------------------------------------------------------------------
            /**
             * evaluates the residual
             * @param[ in ] aResidual matrix to fill with residual
             */
            virtual void compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
            {
                MORIS_ERROR( false, "IWG::compute_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * set the residual size
             * @param[ in ] aResidual matrix to fill with residual
             */
            void set_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
            {
                // FIXME should be 1 if no slave
                aResidual.resize( 2 );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the Jacobian
             * @param[ in ] aJacobians cell of matrices to fill with Jacobians
             */
            virtual void compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
            {
                MORIS_ERROR( false, "IWG::compute_jacobian - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * set the Jacobian
             * @param[ in ] aJacobians cell of matrices to fill with Jacobians
             */
            void set_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
            {
                // get number of dof types for the IWG
                uint tMasterNumDofType = this->get_global_dof_type_list().size();
                uint tSlaveNumDofType  = this->get_global_dof_type_list( mtk::Master_Slave::SLAVE ).size();
                uint tNumDofType = tMasterNumDofType + tSlaveNumDofType;

                // FIXME should be 1 if no slave
                aJacobians.resize( 2 );
                aJacobians( 0 ).resize( tNumDofType );
                aJacobians( 1 ).resize( tNumDofType );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the residual and the Jacobian
             * @param[ in ] aResidual matrix to fill with residual
             * @param[ in ] aJacobians cell of matrices to fill with Jacobians
             */
            virtual void compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                        moris::Cell< Matrix< DDRMat > >                & aResidual )
            {
                MORIS_ERROR( false, " IWG::compute_jacobian_and_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the Jacobian by finite difference
             * @param[ in ] aJacobiansFD  cell of cell of matrices to fill with Jacobians evaluated by FD
             * @param[ in ] aPerturbation real to perturb for FD
             * @param[ in ] aIsMaster     enum master or slave
             */
            virtual void compute_jacobian_FD( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobiansFD,
                                              real                                             aPerturbation,
                                              mtk::Master_Slave                                aIsMaster = mtk::Master_Slave::MASTER )
            {
                // get master and slave number of dof types
                uint tMasterNumDofType = mMasterGlobalDofTypes.size();
                uint tSlaveNumDofType  = mSlaveGlobalDofTypes.size();
                uint tNumDofType       = tMasterNumDofType + tSlaveNumDofType;

                // set the jacobian size
                if ( tSlaveNumDofType > 0 )
                {
                    aJacobiansFD.resize( 2 );
                    aJacobiansFD( 0 ).resize( tNumDofType );
                    aJacobiansFD( 1 ).resize( tNumDofType );
                }
                else
                {
                    aJacobiansFD.resize( 1 );
                    aJacobiansFD( 0 ).resize( tNumDofType );
                }

                // get the number of master residual dofs
                uint tMasterResNumDof = this->get_field_interpolators()( 0 )->get_number_of_space_time_coefficients();
                uint tSlaveResNumDof  = this->get_field_interpolators( mtk::Master_Slave::SLAVE )( 0 )->get_number_of_space_time_coefficients();

                // loop over the IWG dof types
                for( uint iFI = 0; iFI < tNumDofType; iFI++ )
                {
                    // number of dofs for dof type wrt which derivative is computed
                    uint tDerNumDof;

                    // field interpolator wrt which derivative is computed
                    Field_Interpolator* tFI;

                    // coefficients for dof type wrt which derivative is computed
                    Matrix< DDRMat > tCoeff;

                    // derivative wrt master dof type
                    if ( iFI < tMasterNumDofType )
                    {
                        // get field interpolator for derivative
                        tFI = this->get_field_interpolators()( iFI );

                        // get number of master dofs for derivative
                        tDerNumDof = tFI->get_number_of_space_time_coefficients();

                        // get coefficients for derivative
                        tCoeff = tFI->get_coeff();

                        // set size of the Jacobian
                        aJacobiansFD( 0 )( iFI ).set_size( tMasterResNumDof, tDerNumDof, 0.0 );
                        if( tSlaveNumDofType > 0 )
                        {
                            aJacobiansFD( 1 )( iFI ).set_size( tSlaveResNumDof, tDerNumDof, 0.0 );
                        }
                    }
                    // derivative wrt slave dof type
                    else
                    {
                        // get field interpolator for derivative
                        tFI = this->get_field_interpolators( mtk::Master_Slave::SLAVE )( iFI - tMasterNumDofType );

                        // get number of slave dofs for derivative
                        tDerNumDof = tFI->get_number_of_space_time_coefficients();

                        // get coefficients for derivative
                        tCoeff = tFI->get_coeff();

                        // set size of the Jacobian
                        aJacobiansFD( 0 )( iFI ).set_size( tMasterResNumDof, tDerNumDof, 0.0 );
                        aJacobiansFD( 1 )( iFI ).set_size( tSlaveResNumDof,  tDerNumDof, 0.0 );
                    }

                    for( uint iCoeff = 0; iCoeff < tDerNumDof; iCoeff++ )
                    {
                        // perturbation of the coefficent
                        Matrix< DDRMat > tCoeffPert = tCoeff;
                        tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation;

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // evaluate the residual
                        moris::Cell< Matrix< DDRMat > > tResidual_Plus;
                        this->compute_residual( tResidual_Plus );

                        // perturbation of the coefficent
                        tCoeffPert = tCoeff;
                        tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation;

                        // setting the perturbed coefficients
                        tFI->set_coeff( tCoeffPert );

                        // evaluate the residual
                        moris::Cell< Matrix< DDRMat > > tResidual_Minus;
                        this->compute_residual( tResidual_Minus );

                        // evaluate Jacobian
                        aJacobiansFD( 0 )( iFI ).get_column( iCoeff ) = ( tResidual_Plus( 0 ) - tResidual_Minus( 0 ) )/ ( 2.0 * aPerturbation );
                        if( tSlaveNumDofType > 0 )
                        {
                            aJacobiansFD( 1 )( iFI ).get_column( iCoeff ) = ( tResidual_Plus( 1 ) - tResidual_Minus( 1 ) )/ ( 2.0 * aPerturbation );
                        }
                    }
                    // reset the coefficients values
                    tFI->set_coeff( tCoeff );
                }
            }

//------------------------------------------------------------------------------
//
//            virtual real compute_integration_error( const Matrix< DDRMat > & aNodalDOF,
//                                                    real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
//                                                    const uint        & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }

//------------------------------------------------------------------------------
//
//            real interpolate_scalar_at_point( const Matrix< DDRMat > & aNodalWeakBC,
//                                              const uint             & aPointIndex )
//            {
//                MORIS_ERROR( false, "This function does nothing" );
//                return 0.0;
//            }


        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_IWG_HPP_ */
