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

            // active dof types
            moris::Cell< moris::Cell< MSI::Dof_Type > > mActiveDofTypes;

            // field interpolators
            moris::Cell< Field_Interpolator* > mMasterFI;
            moris::Cell< Field_Interpolator* > mSlaveFI;

            // active property type
            moris::Cell< fem::Property_Type > mActivePropertyTypes;

            // properties
            moris::Cell< Property* > mMasterProp;

            // global active dof types
            moris::Cell< moris::Cell< MSI::Dof_Type > > mActiveDofTypesGlobal;

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
             * set field interpolators
             * @param[ in ] aFieldInterpolators cell of field interpolator pointers
             */
            void set_field_interpolators( moris::Cell< Field_Interpolator* > & aFieldInterpolators,
                                          mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                // check input size
                MORIS_ASSERT( aFieldInterpolators.size() == mActiveDofTypes.size(), "IWG::set_field_interpolators - wrong input size. " );

                // set field interpolators
                if ( aIsMaster == mtk::Master_Slave::MASTER )
                {
                    // set master field interpolators
                    mMasterFI = aFieldInterpolators;
                }
                else
                {
                    // set slave field interpolators
                    mSlaveFI = aFieldInterpolators;
                }

            }

            /**
              * check that field interpolators were assigned
              */
             void check_field_interpolators( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // master field interpolators
                 if ( aIsMaster == mtk::Master_Slave::MASTER )
                 {
                     // check master field interpolators cell size
                     MORIS_ASSERT( mMasterFI.size() == mActiveDofTypes.size(), "IWG::check_field_interpolators - wrong master FI size. " );

                     // loop over all master FI and check that they are assigned
                     for( uint iFI = 0; iFI < mActiveDofTypes.size(); iFI++ )
                     {
                         MORIS_ASSERT( mMasterFI( iFI ) != nullptr, "IWG::check_field_interpolators - master FI missing. " );
                     }
                 }
                 // slave field interpolators
                 else
                 {
                    // check slave field interpolators cell size
                    MORIS_ASSERT( mSlaveFI.size() == mActiveDofTypes.size(), "IWG::check_field_interpolators - wrong slave FI size. " );

                    // loop over all slave FI and check that they are assigned
                    for( uint iFI = 0; iFI < mActiveDofTypes.size(); iFI++ )
                    {
                        MORIS_ASSERT( mSlaveFI( iFI ) != nullptr, "IWG::check_field_interpolators - slave FI missing. " );
                    }
                 }
             }

//------------------------------------------------------------------------------
             /**
              * set properties
              * @param[ in ] aProperties cell of property pointers
              */
             void set_properties( moris::Cell< Property* > & aProperties,
                                  mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
             {
                 // check input size
                 MORIS_ASSERT( aProperties.size() == mActivePropertyTypes.size(), "IWG::set_properties - wrong input size. " );

                 // fixme check we pass in the right properties

                 // set properties
                 mMasterProp = aProperties;

                 // create a global dof type list
                 this->create_global_dof_types_list();
             }

             /**
               * check that properties were assigned
               */
              void check_properties( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
              {
                  // master properties
                  if ( aIsMaster == mtk::Master_Slave::MASTER )
                  {
                      // check master field interpolators cell size
                      MORIS_ASSERT( mMasterProp.size() == mActivePropertyTypes.size(), "IWG::check_properties - wrong master Prop size. " );

                      // loop over all master FI and check that they are assigned
                      for( uint iProp = 0; iProp < mActivePropertyTypes.size(); iProp++ )
                      {
                          MORIS_ASSERT( mMasterProp( iProp ) != nullptr, "IWG::check_properties - master property missing. " );
                      }
                  }
              }

//------------------------------------------------------------------------------

              void create_global_dof_types_list()
              {
                  // set the size of the active dof type list
                  uint tCounterMax = this->get_active_dof_types().size();

                  for ( Property* tProperty : mMasterProp )
                  {
                      tCounterMax += tProperty->get_active_dof_types().size();
                  }
                  mActiveDofTypesGlobal.resize( tCounterMax );
                  moris::Cell< sint > tCheckList( tCounterMax, -1 );

                  // init total dof counter
                  uint tCounter = 0;

                  // get active dof type for IWG
                  for ( uint iDOF = 0; iDOF < mActiveDofTypes.size(); iDOF++ )
                  {
                      tCheckList( tCounter ) = static_cast< uint >( mActiveDofTypes( iDOF )( 0 ) );
                      mActiveDofTypesGlobal( tCounter ) = mActiveDofTypes( iDOF );
                      tCounter++;
                  }

                  for ( Property* tProperty : mMasterProp )
                  {
                      // get active dof type
                      moris::Cell< moris::Cell< MSI::Dof_Type > > tActiveDofType = tProperty->get_active_dof_types();

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
                              tCheckList( tCounter ) = static_cast< uint >( tActiveDofType( iDOF )( 0 ) );
                              mActiveDofTypesGlobal( tCounter ) = tActiveDofType( iDOF );
                              tCounter++;
                          }
                      }
                  }

                  // get the number of unique dof type groups, i.e. the number of interpolators
                  mActiveDofTypesGlobal.resize( tCounter );
              };

              const moris::Cell< moris::Cell< MSI::Dof_Type > > & get_global_dof_types_list() const
              {
                  return mActiveDofTypesGlobal;
              };

//------------------------------------------------------------------------------
            /**
             * returns a dof type for the residual
             */
            moris::Cell< MSI::Dof_Type > get_residual_dof_type() const
            {
                return mResidualDofType;
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of dof types used to evaluate the residual
             * and the jacobian
             */
            moris::Cell< moris::Cell< MSI::Dof_Type > > get_active_dof_types() const
            {
                return mActiveDofTypes;
            };

//------------------------------------------------------------------------------
            /**
             * returns a cell of mp types used to evaluate the residual
             * and the jacobian
             */
            moris::Cell< fem::Property_Type > get_active_property_types() const
            {
                return mActivePropertyTypes;
            };

//------------------------------------------------------------------------------
            /**
             * evaluates the residual
             */
            virtual void compute_residual( Matrix< DDRMat > & aResidual )
            {
                MORIS_ERROR( false, "IWG::compute_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the Jacobian
             */
            virtual void compute_jacobian( moris::Cell< Matrix< DDRMat > > & aJacobians )
            {
                MORIS_ERROR( false, "IWG::compute_jacobian - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the residual and the Jacobian
             */
            virtual void compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > > & aJacobians,
                                                        Matrix< DDRMat >                & aResidual )
            {
                MORIS_ERROR( false, " IWG::compute_jacobian_and_residual - This function does nothing. " );
            }

//------------------------------------------------------------------------------
            /**
             * evaluates the Jacobian by finite difference
             */
            virtual void compute_jacobian_FD( moris::Cell< Matrix< DDRMat > > & aJacobiansFD,
                                              real                              aPerturbation,
                                              mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
            {
                if( aIsMaster == mtk::Master_Slave::MASTER )
                {
                    // set the jacobian size
                    aJacobiansFD.resize( mMasterFI.size() );

                    // loop over the left field interpolator
                    for( uint iFI = 0; iFI < mMasterFI.size(); iFI++ )
                    {
                        aJacobiansFD( iFI ).set_size( mMasterFI( 0 )->get_number_of_space_time_coefficients(),
                                                      mMasterFI( iFI )->get_number_of_space_time_coefficients(),
                                                      0.0 );

                        // get field interpolator coefficients
                        Matrix< DDRMat > tCoeff = mMasterFI( iFI )->get_coeff();

                        for( uint iCoeff = 0; iCoeff< mMasterFI( iFI )->get_number_of_space_time_coefficients(); iCoeff++ )
                        {
                            // perturbation of the coefficent
                            Matrix< DDRMat > tCoeffPert = tCoeff;
                            tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation;

                            // setting the perturbed coefficients
                            mMasterFI( iFI )->set_coeff( tCoeffPert );

                            // evaluate the residual
                            Matrix< DDRMat > tResidual_Plus;
                            this->compute_residual( tResidual_Plus);

                            // perturbation of the coefficent
                            tCoeffPert = tCoeff;
                            tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation;

                            // setting the perturbed coefficients
                            mMasterFI( iFI )->set_coeff( tCoeffPert );

                            // evaluate the residual
                            Matrix< DDRMat > tResidual_Minus;
                            this->compute_residual( tResidual_Minus );

                            // evaluate Jacobian
                            aJacobiansFD( iFI ).get_column( iCoeff ) = ( tResidual_Plus - tResidual_Minus )/ ( 2.0 * aPerturbation );

                        }
                        // reset the coefficients values
                        mMasterFI( iFI )->set_coeff( tCoeff );
                    }
                }
                else
                {
                    // set the jacobian size
                    aJacobiansFD.resize( mSlaveFI.size() );

                    // loop over the right field interpolators
                    for( uint iRightFI = 0; iRightFI < mSlaveFI.size(); iRightFI++ )
                    {
                        // get field interpolator coefficients
                        Matrix< DDRMat > tCoeff = mSlaveFI( iRightFI )->get_coeff();

                        for( uint iCoeff = 0; iCoeff< mSlaveFI( iRightFI )->get_number_of_space_time_coefficients(); iCoeff++ )
                        {
                            // perturbation of the coefficent
                            Matrix< DDRMat > tCoeffPert = tCoeff;
                            tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) + aPerturbation;

                            // setting the perturbed coefficients
                            mSlaveFI( iRightFI )->set_coeff( tCoeffPert );

                            // evaluate the residual
                            Matrix< DDRMat > tResidual_Plus;
                            this->compute_residual( tResidual_Plus );

                            // perturbation of the coefficent
                            tCoeffPert = tCoeff;
                            tCoeffPert( iCoeff ) = tCoeffPert( iCoeff ) - aPerturbation;

                            // setting the perturbed coefficients
                            mSlaveFI( iRightFI )->set_coeff( tCoeffPert );

                            // evaluate the residual
                            Matrix< DDRMat > tResidual_Minus;
                            this->compute_residual( tResidual_Minus );

                            // evaluate Jacobian
                            aJacobiansFD( iRightFI ).get_column( mMasterFI( iRightFI )->get_number_of_space_time_coefficients() + iCoeff )
                                = ( tResidual_Plus - tResidual_Minus )/ ( 2.0 * aPerturbation );

                        }
                        // reset the coefficients values
                        mSlaveFI( iRightFI )->set_coeff( tCoeff );
                    }
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
