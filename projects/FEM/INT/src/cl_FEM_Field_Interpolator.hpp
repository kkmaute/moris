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

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_MSI_Dof_Type_Enums.hpp"     //FEM/MSI/src

#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
class Property;

        class Field_Interpolator
        {
            // how many fields are to be interpolated
            const uint mNumberOfFields;

            // pointer to space and time interpolation objects
            Interpolation_Function_Base * mSpaceInterpolation = nullptr;
            Interpolation_Function_Base * mTimeInterpolation  = nullptr;

            // space and time geometry interpolator
            Geometry_Interpolator * mGeometryInterpolator = nullptr;

            // space, time, and space time number of bases
            uint mNSpaceBases;
            uint mNTimeBases;
            uint mNFieldBases;

            // space time number of coefficients
            uint mNFieldCoeff;

            // space, time dimensions
            uint mNSpaceDim;
            uint mNTimeDim;

            // space parametric dimensions
            uint mNSpaceParamDim;

            // parametric point where field is interpolated
            Matrix< DDRMat > mXi;
            Matrix< DDRMat > mTau;

            // matrix of field coefficients uHat
            Matrix < DDRMat > mUHat;

            // dof type
            moris::Cell< MSI::Dof_Type > mDofType;

            // dv type
            moris::Cell< PDV_Type > mDvType;

            // flag for evaluation
            bool mNBuildEval = true;
            bool mNEval      = true;
            bool mdNdxEval   = true;
            bool md2Ndx2Eval = true;
            bool md3Ndx3Eval = true;
            bool mdNdtEval   = true;
            bool md2Ndt2Eval = true;
            bool md2NdxtEval = true;

            // storage
            Matrix< DDRMat > mNBuild;
            Matrix< SDRMat > mN;
            Matrix< DDRMat > mdNdx;
            Matrix< DDRMat > md2Ndx2;
            Matrix< DDRMat > md3Ndx3;
            Matrix< DDRMat > mdNdt;
            Matrix< DDRMat > md2Ndt2;
            Matrix< DDRMat > md2Ndxt;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aNumberOfFields           number of interpolated fields
             * @param[ in ] aFieldInterpolationRule   field interpolation rule
             * @param[ in ] aGeometryInterpolator     pointer to geometry interpolator object
             * @param[ in ] aDofType                  dof type for the interpolated fields
             */
            Field_Interpolator( const uint                         & aNumberOfFields,
                                const Interpolation_Rule           & aFieldInterpolationRule,
                                      Geometry_Interpolator*         aGeometryInterpolator,
                                const moris::Cell< MSI::Dof_Type >   aDofType );

            /**
             * constructor
             * @param[ in ] aNumberOfFields           number of interpolated fields
             * @param[ in ] aFieldInterpolationRule   field interpolation rule
             * @param[ in ] aGeometryInterpolator     pointer to geometry interpolator object
             * @param[ in ] aDvType                   dv type for the interpolated fields
             */
            Field_Interpolator( const uint                         & aNumberOfFields,
                                const Interpolation_Rule           & aFieldInterpolationRule,
                                      Geometry_Interpolator*         aGeometryInterpolator,
                                const moris::Cell< PDV_Type >          aDvType );

            /**
             * trivial constructor for unit test
             */
            Field_Interpolator( const uint & aNumberOfFields,
                                const moris::Cell< MSI::Dof_Type >   aDofType ) : mNumberOfFields( aNumberOfFields ),
                                                                                  mDofType( aDofType )
            {
                mNFieldCoeff = mNumberOfFields;
            };

            /**
              * trivial constructor for unit test
              */
             Field_Interpolator( const uint                  & aNumberOfFields,
                                 const moris::Cell< PDV_Type >   aDvType ) : mNumberOfFields( aNumberOfFields ),
                                                                                 mDvType( aDvType )
             {
                 mNFieldCoeff = mNumberOfFields;
             };

//------------------------------------------------------------------------------
            /**
             * default constructor
             */
            ~Field_Interpolator();

            //------------------------------------------------------------------------------
            /**
             * reset evaluation flags
             */
            void reset_eval_flags();

//------------------------------------------------------------------------------
            /**
             * get dof type
             */
            const moris::Cell< MSI::Dof_Type > & get_dof_type() const
            {
                return mDofType;
            }

//------------------------------------------------------------------------------
            /**
             * get dof type
             */
            const moris::Cell< PDV_Type > & get_dv_type() const
            {
                return mDvType;
            }

//------------------------------------------------------------------------------
            /**
             * get the number of fields
             */
            uint get_number_of_fields() const
            {
                return mNumberOfFields;
            }

//------------------------------------------------------------------------------
            /**
             * get the number of space bases
             */
            uint get_number_of_space_bases() const
            {
                return mNSpaceBases;
            }

//------------------------------------------------------------------------------
            /**
             * get the number of time bases
             */
            uint get_number_of_time_bases() const
            {
                return mNTimeBases;
            }

//------------------------------------------------------------------------------
             /**
              * get the number of space time bases
              */
             uint get_number_of_space_time_bases() const
             {
                 return mNFieldBases;
             }

//------------------------------------------------------------------------------
             /**
              * get the number of space time coefficients
              */
             uint get_number_of_space_time_coefficients() const
             {
                 return mNFieldCoeff;
             }
//------------------------------------------------------------------------------
             /**
              * get space interpolation order
              */
             mtk::Interpolation_Order get_space_interpolation_order()
             {
                 return mSpaceInterpolation->get_interpolation_order();
             }
//------------------------------------------------------------------------------
            /**
             * set the parametric point where field is interpolated
             * @param[ in ] aParamPoint evaluation point in space and time
             */
            void set_space_time( const Matrix< DDRMat > & aParamPoint );

//------------------------------------------------------------------------------
             /**
              * get the parametric point in space where field is interpolated
              */
             const Matrix< DDRMat > & get_space() const
             {
                 return mXi;
             }

//------------------------------------------------------------------------------
            /**
             * get the parametric point in time where field is interpolated of tau
             */
            const Matrix< DDRMat > & get_time() const
            {
                return mTau;
            }

//------------------------------------------------------------------------------
             /**
              * set the coefficients of the field uHat
              * @param[ in ] aUHat coefficients
              */
             void set_coeff( const Matrix< DDRMat > & aUHat );

//------------------------------------------------------------------------------
             /**
             * get the coefficients of the field uHat
             * @param[ out ] mUHat coefficients
             */
              const Matrix< DDRMat > & get_coeff() const
              {
                  return mUHat;
              }

//------------------------------------------------------------------------------
            /**
             * return the space time shape functions
             * @param[ out ] shape functions matrix
             *               ( 1 x <number of bases> )
             */
            const Matrix < DDRMat > & NBuild();

//------------------------------------------------------------------------------
            /**
             * evaluate the space time shape functions
             */
            void eval_NBuild();

//------------------------------------------------------------------------------
            /**
             * return the N for vectorial field ( space time shape functions )
             * @param[ out ] ( nNumberOfFields x mNFieldCoeff )
             */
             const Matrix < SDRMat > & N();

//------------------------------------------------------------------------------
            /**
             * evaluate the N for vectorial field ( space time shape functions )
             * @param[ out ] ( nNumberOfFields x mNFieldCoeff )
             */
             void eval_N();

//------------------------------------------------------------------------------
            /**
             * return the nth order derivatives of the space time shape functions
             * wrt space x
             * @param[ in ]  aDerivativeOrder derivative order
             * @param[ out ] dnNdxn           nth order spatial derivative of the shape functions
             */
            const Matrix< DDRMat > & dnNdxn( const uint & aDerivativeOrder );

//------------------------------------------------------------------------------
            /**
             * evaluates the first derivatives of the space time shape functions
             * wrt space x
             */
            void eval_d1Ndx1();

//------------------------------------------------------------------------------
            /**
             * evaluates the second derivatives of the space time shape functions
             * wrt space x
             */
            void eval_d2Ndx2();

//------------------------------------------------------------------------------
            /**
             * evaluates the third derivatives of the space time shape functions
             * wrt space x
             */
            void eval_d3Ndx3();

//------------------------------------------------------------------------------
            /**
             * return the nth order derivatives of the space time shape functions
             * wrt time t
             * @param[ in ]  aDerivativeOrder derivative order
             * @param[ out ] dnNdtn           nth order time derivative of the shape functions
             */
            const Matrix< DDRMat > & dnNdtn( const uint & aDerivativeOrder );

//------------------------------------------------------------------------------
            /**
             * evaluates the first derivative of the space time shape functions
             * wrt time t
             */
            void eval_d1Ndt1();

//------------------------------------------------------------------------------
            /**
            * evaluates the second derivative of the space time shape functions
            * wrt time t
            */
            void eval_d2Ndt2();

//------------------------------------------------------------------------------
            /**
             * return the 1st order mixed derivatives of the space time shape functions
             * @param[ out ] d2Ndxt  mixed space & time derivative of the shape functions
             */
            const Matrix< DDRMat > & d2Ndxt();

//------------------------------------------------------------------------------
            /**
             * mixed space & time derivative of the space time shape functions
             * wrt time x & t
             */
            void eval_d2Ndxt();

//------------------------------------------------------------------------------
            /**
            * evaluates the field at given space and time Xi, Tau
            * @param[ out ]          interpolated field
            */
            Matrix< DDRMat > val();

//------------------------------------------------------------------------------
            /**
            * evaluates the field space derivatives at given space and time evaluation point
            * @param[ in ]  aDerivativeOrder  order of the required derivatives
            * @param[ out ] gradx             space derivatives
            */
            Matrix< DDRMat > gradx( const uint & aDerivativeOrder );

//------------------------------------------------------------------------------
            /**
            * evaluates the field spatial divergence at given space and time evaluation point
            * @return divergence of the field
            */
            moris::real div();

            /**
            * evaluates the spatial divergence operator at given space and time evaluation point
            * @param[ out ] divergence operator
            */
            Matrix< DDRMat > div_operator();

//------------------------------------------------------------------------------
            /**
             * evaluates the field time derivative at given space and time evaluation point
             * @param[ in ] aDerivativeOrder  order of the required derivative
             * @param[ out ] gradt            time derivatives
             */
            Matrix< DDRMat > gradt( const uint & aDerivativeOrder );

//------------------------------------------------------------------------------
            /**
             * evaluates the mixed field space & time derivative at given space and time evaluation point
             * @param[ out ] gradxt            mixed space time derivatives
             */
            Matrix< DDRMat > gradxt();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_ */
