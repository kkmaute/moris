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

            // active field interpolators
        	//moris::Cell< MSI::Dof_Type > mDofDependencies;
            //moris::Cell< Field_Interpolator* > mDependencyFieldInterpolators;

            // a function to set the model parameter coefficients
            real ( *mFunc )( Matrix< DDRMat > aPhysPoint );


//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             */
//            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
//                                          const Interpolation_Rule     & aFieldInterpolationRule,
//                                          const Geometry_Interpolator*   aGeometryInterpolator,
//                                          fem::Mp_Type                   aMpType )
//                                        : mGeometryInterpolator( aGeometryInterpolator ){};

//            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
//                                          const Interpolation_Rule     & aFieldInterpolationRule,
//                                          const Geometry_Interpolator*   aGeometryInterpolator,
//                                          MSI::Dof_Type                  aDofType )
//                                        : mGeometryInterpolator( aGeometryInterpolator ){};

//            Model_Parameter_Interpolator( const uint                    aNumberOfFields,
//                                          const Interpolation_Rule     & aFieldInterpolationRule,
//                                          const Geometry_Interpolator*   aGeometryInterpolator,
//                                          MSI::Pdv_Type                  aPdvType )
//                                        : mGeometryInterpolator( aGeometryInterpolator ){};

            Model_Parameter_Interpolator( const uint                   aNumberOfFields,
                                          const Interpolation_Rule&    aFieldInterpolationRule,
                                          const Geometry_Interpolator* aGeometryInterpolator,
                                          real ( *aFunc)( Matrix< DDRMat > aPhysPoint ) )
            : mGeometryInterpolator( aGeometryInterpolator ),
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

                 // evaluate the coefficients from the provided function
                 // from geometry interpolator of the interpolation cell
                 Matrix< DDRMat > tPhysCoords = mGeometryInterpolator->build_space_time_phys_coords();

                 for( uint iCoeff = 0; iCoeff < tNumCoeffs; iCoeff++ )
                 {
                    tUHat( iCoeff ) = mFunc( tPhysCoords.get_column( iCoeff ) );
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
            * evaluates the derivative of the model parameter wrt its coefficients
            * @param[ out ] first derivative r
            */
           Matrix < DDRMat > N( )
           {
               return mFieldInterpolator->N();
           }

//------------------------------------------------------------------------------
//           /**
//            * evaluates the derivative of the model parameter wrt a given dof type
//            * @param[ out ] first derivative of the model parameter wrt to a dof_type
//            */
//           Matrix < DDRMat > gradDof( const msi::Dof_Type aDofType )
//           {
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
