/*
 * cl_FEM_Field_Interpolator.hpp
 *
 *  Created on: Jan 31, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp" //LNA/src
#include "cl_Cell.hpp"
#include "linalg_typedefs.hpp" //LNA/src
#include "assert.h"

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Field_Interpolator
        {
            // how many fields are to be interpolated
            const uint mNumberOfFields;

//            // boolean set to true if space interpolation only
//            bool mSpaceOnlyFlag;

            // pointer to space and time interpolation objects
            Interpolation_Function_Base * mSpaceInterpolation = nullptr;
            Interpolation_Function_Base * mTimeInterpolation  = nullptr;

            // space and time geometry interpolator objects
            const Geometry_Interpolator * mGeometryInterpolator = nullptr;

            // space, time, and space time number of bases
            uint mNSpaceBases;
            uint mNTimeBases;
            uint mNFieldBases;

            // space time number of coefficients
            uint mNFieldCoeff;

            // space, time dimensions
            uint mNSpaceDim;
            uint mNTimeDim;

            // parametric point where field is interpolated
            Matrix< DDRMat > mXi;
            Matrix< DDRMat > mTau;

            // matrix of field coefficients uHat
            Matrix < DDRMat > mUHat;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aNumberOfFields           number of interpolated fields
             * @param[ in ] aFieldInterpolationRule   pointer to field interpolation rule
             * @param[ in ] aGeometryInterpolator     pointer to geometry interpolator object
             *
             */
            Field_Interpolator( const uint                   & aNumberOfFields,
                                const Interpolation_Rule     & aFieldInterpolationRule,
                                const Geometry_Interpolator*   aGeometryInterpolator );

//------------------------------------------------------------------------------
            /**
             * default constructor
             */
            ~Field_Interpolator();

//------------------------------------------------------------------------------
            /**
             * get the number of fields
             */
            uint const get_number_of_fields() const
            {
                return mNumberOfFields;
            }

//------------------------------------------------------------------------------
            /**
             * get the number of space bases
             */
            uint const get_number_of_space_bases() const
            {
                return mNSpaceBases;
            }

//------------------------------------------------------------------------------
            /**
             * get the number of time bases
             */
            uint const get_number_of_time_bases() const
            {
                return mNTimeBases;
            }

//------------------------------------------------------------------------------
             /**
              * get the number of space time bases
              */
             uint const get_number_of_space_time_bases() const
             {
                 return mNFieldBases;
             }

//------------------------------------------------------------------------------
             /**
              * get the number of space time coefficients
              */
             uint const get_number_of_space_time_coefficients() const
             {
                 return mNFieldCoeff;
             }
//------------------------------------------------------------------------------
            /**
             * set the parametric point where field is interpolated xi, tau
             */
            void set_space_time( const Matrix< DDRMat > & aParamPoint );

//------------------------------------------------------------------------------
            /**
             * set the parametric point where field is interpolated xi
             */
            void set_space( const Matrix< DDRMat > & aXi );

//------------------------------------------------------------------------------
             /**
              * get the parametric point where field is interpolated of xi, tau
              */
              Matrix< DDRMat > get_space() const
              {
                  return mXi;
              }

              Matrix< DDRMat > get_time() const
              {
                  return mTau;
              }

//------------------------------------------------------------------------------
             /**
              * set the coefficients of the field uHat
              */
             void set_coeff( const Matrix< DDRMat > & aUHat );

//------------------------------------------------------------------------------
             /**
             * get the coefficients of the field uHat
             */
              Matrix< DDRMat > get_coeff() const
              {
                  return mUHat;
              }

////------------------------------------------------------------------------------
//             /**
//              * get the space only flag
//              */
//              const bool space_only() const
//              {
//                   return mSpaceOnlyFlag;
//              }

//------------------------------------------------------------------------------
            /**
             * evaluates the space time shape functions
             * @param[ out ] shape function matrix
             *             ( 1 x <number of basis> )
             */
            Matrix < DDRMat > N();

//------------------------------------------------------------------------------
            /**
             * evaluates the first derivatives of the space time shape functions
             * wrt space x
             * @param[ out ] dNdx
             *             ( < number of space dimensions > x <number of space time basis > )
             */
             Matrix< DDRMat > Bx();

//------------------------------------------------------------------------------
            /**
             * evaluates the second derivatives of the space time shape functions
             * wrt space x
             * @param[ out ]         d2Ndx2
             *             ( < number of space dimensions > x <number of space time basis > )
             */
            Matrix< DDRMat > eval_d2Ndx2();

//------------------------------------------------------------------------------
            /**
             * evaluates the first derivative of the space time shape functions
             * wrt time t
             * @param[ out ]       dNdt
             *             ( < number of time dimensions > x <number of space time basis > )
             */
            Matrix< DDRMat > Bt();

//------------------------------------------------------------------------------
            /**
            * evaluates the second derivative of the space time shape functions
            * wrt time t
            * @param[ out ]       d2Ndt2
            *             ( < number of time dimensions > x <number of space time basis > )
            */
            Matrix< DDRMat > eval_d2Ndt2();

//------------------------------------------------------------------------------
            /**
            * evaluates the field at given space and time Xi, Tau
            * @param[ out ]          interpolated field
            */
            Matrix< DDRMat > val();

//------------------------------------------------------------------------------
            /**
            * evaluates the field space derivative at given space and time Xi, Tau
            * @param[ in ] aDerivativeOrder  order of the required derivative
            * @param[ in ] aXHat             space coordinates of the treated element
            *             ( < number of space dimensions > x < number of space basis > )
            *
            */
            Matrix< DDRMat > gradx( const uint & aDerivativeOrder );

//------------------------------------------------------------------------------
            /**
             * evaluates the field time derivative at given space and time Xi, Tau
             * @param[ in ] aDerivativeOrder  order of the required derivative
             * @param[ in ] aTHat  time coordinates of the treated element
             *             ( < number of time dimensions > x < number of time basis > )
             */
            Matrix< DDRMat > gradt( const uint & aDerivativeOrder );

//------------------------------------------------------------------------------
            /**
             * evaluates the determinant of the Jacobian mapping
             * at given space and time Xi, Tau
             */
            real det_J();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_ */
