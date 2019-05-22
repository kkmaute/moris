/*
 * cl_FEM_Model_Parameter_Interpolator.hpp
 *
 *  Created on: May 02, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_MODEL_PARAMETER_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_MODEL_PARAMETER_INTERPOLATOR_HPP_

#include "typedefs.hpp"                  //MRS/COR/src

#include "cl_Matrix.hpp"                 //LNA/src
#include "cl_Cell.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_Field_Interpolator.hpp" //FEM/src
#include "cl_MSI_Dof_Type_Enums.hpp" //FEM/MSI/src


namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Model_Parameter_Interpolator
        {

        private:

            // field interpolator pointer associated with the model parameter
            Field_Interpolator* mFieldInterpolator = nullptr;

            // geometry interpolator associated with the model parameter
            const Geometry_Interpolator* mGeometryInterpolator = nullptr;

            // field interpolator type
            fem::Mp_Type mMPType;

            // active field interpolators
            moris::Cell< MSI::Dof_Type > mDofDependencies;

            //moris::Cell< Field_Interpolator* > mDependencyFieldInterpolators;

            // a function to set the model parameter coefficients
            real ( *mFunc )( Matrix< DDRMat > aSpacePhysPoint,
                             Matrix< DDRMat > aTimePhysPoint );


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             */
//            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
//                                          const Interpolation_Rule     & aFieldInterpolationRule,
//                                          const Geometry_Interpolator*   aGeometryInterpolator,
//                                          fem::MP_Type                   aMpType )
//                                        : mGeometryInterpolator( aGeometryInterpolator ){};

//            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
//                                          const Interpolation_Rule    & aFieldInterpolationRule,
//                                          const Geometry_Interpolator * aGeometryInterpolator,
//                                          MSI::Dof_Type                 aDofType )
//                                        : mGeometryInterpolator( aGeometryInterpolator ){};

//            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
//                                          const Interpolation_Rule     & aFieldInterpolationRule,
//                                          const Geometry_Interpolator*   aGeometryInterpolator,
//                                          MSI::Pdv_Type                  aPdvType )
//                                        : mGeometryInterpolator( aGeometryInterpolator ){};

            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
                                          const Interpolation_Rule    & aFieldInterpolationRule,
                                          const Geometry_Interpolator * aGeometryInterpolator,
                                          fem::Mp_Type                  aMpType,
                                          real ( *aFunc)( Matrix< DDRMat > aSpacePhysPoint,
                                                          Matrix< DDRMat > aTimePhysPoint ) )
            : mGeometryInterpolator( aGeometryInterpolator ),
              mMPType( aMpType ),
              mFunc( aFunc )
            {
                // create a field interpolator pointer
                mFieldInterpolator = new Field_Interpolator( aNumberOfFields,
                                                             aFieldInterpolationRule,
                                                             mGeometryInterpolator );
            };

//------------------------------------------------------------------------------
            /**
             * default destructor
             */
            ~Model_Parameter_Interpolator()
            {
                if ( mFieldInterpolator != NULL )
                {
                    delete mFieldInterpolator;
                }
            }

//------------------------------------------------------------------------------
            /**
             * set the parametric point where model parameter is interpolated
             * @param[ in ] aParamPoint evaluation point in space and time
             */
            void set_space_time( const Matrix< DDRMat > & aParamPoint )
            {
                mFieldInterpolator->set_space_time( aParamPoint );
            };

//------------------------------------------------------------------------------
             /**
              * set the coefficients of the model parameter u
              * @param[ in ] aUHat coefficients
              */
             void set_coeff( const Matrix< DDRMat > & aUHat )
             {
                 mFieldInterpolator->set_coeff( aUHat );
             }

             void set_coeff()
             {
                 // init the coefficient vector
                 uint tNumCoeffs = mFieldInterpolator->get_number_of_space_time_coefficients();
                 Matrix< DDRMat > tUHat( tNumCoeffs, 1, 0.0 );


                 // from geometry interpolator of the interpolation cell
                 Matrix< DDRMat > tSpacePhysCoords = mGeometryInterpolator->get_space_coeff();
                 Matrix< DDRMat > tTimePhysCoords  = mGeometryInterpolator->get_time_coeff();

                 // evaluate the coefficients from the provided function
                 uint tCoeffCounter = 0;
                 for( uint iTime = 0; iTime < tTimePhysCoords.n_rows(); iTime++ )
                 {
                     for( uint iSpace = 0; iSpace < tSpacePhysCoords.n_rows(); iSpace++ )
                     {
                         tUHat( tCoeffCounter ) = mFunc( tSpacePhysCoords.get_row( iSpace ),
                                                         tTimePhysCoords.get_row( iTime ) );
                         tCoeffCounter++;
                     }
                 }

                 // set the coefficients for the field interpolator
                 mFieldInterpolator->set_coeff( tUHat );
             }

//------------------------------------------------------------------------------
             /**
             * get the coefficients of the model parameter uHat
             * @param[ out ] mUHat coefficients
             */
              Matrix< DDRMat > get_coeff() const
              {
                  return mFieldInterpolator->get_coeff();
              }

//------------------------------------------------------------------------------
            /**
            * evaluates the model parameter at given point
            * @param[ out ] interpolated model parameter
            */
            Matrix< DDRMat > val()
            {
                return mFieldInterpolator->val();
            };

//------------------------------------------------------------------------------

            Matrix< DDRMat > gradx( const uint aDerivativeOrder )
            {
                return mFieldInterpolator->gradx( aDerivativeOrder );
            }

//------------------------------------------------------------------------------

            Matrix< DDRMat > gradt( const uint aDerivativeOrder )
            {
                return mFieldInterpolator->gradt( aDerivativeOrder );
            }

//------------------------------------------------------------------------------
           /**
            * evaluates the derivative of the model parameter wrt a given dof type
            * @param[ out ] first derivative of the model parameter wrt to a dof_type
            */
           Matrix < DDRMat > val_gradDof( MSI::Dof_Type aDofType )
           {
               Matrix< DDRMat > tValGradDof( mFieldInterpolator->get_number_of_fields(),
                                             mFieldInterpolator->get_number_of_space_time_coefficients(),
                                             0.0 );

               // if no dof dependencies
               if ( mDofDependencies.size()>0 )
               {
                   MORIS_ASSERT( false, "Model_Parameter_Interpolator::Modeval_gradDof - not implemented when dependencies.");
               }
               return tValGradDof;
           }

//           /**
//             * evaluates the derivative of the model parameter wrt a given dof type
//             * @param[ out ] first derivative of the model parameter wrt to a dof_type
//             */
//            Matrix < DDRMat > gradDof( const msi::Dof_Type aDofType )
//            {
//        	   // fixme two cases created with a func or with a model parameter function
//
//               // get the derivative of the model parameter wrt to its coeff
//               Matrix< DDRMat > tN = this->N() ;
//
//               // get the potential derivative of the model parameter coefficients
//               // wrt the dof type
//               // cascade down to lower level???
//               Matrix< DDRMat > tGradDependencies;
//               for ( uint iDep = 0; iDep < mDofDependencies.size(); iDep++ )
//               {
//                   tGradDofDependencies = tGradDofDependencies
//                       + mFuncDeriv->val() * mDependencyFieldInterpolators( iDep )->gradDof( aDofType );
//               }
//
//               // return the derivative
//               return tN * tGradDofDependencies;
//           }
//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_MODEL_PARAMETER_BASE_HPP_ */
