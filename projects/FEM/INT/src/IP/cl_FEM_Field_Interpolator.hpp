/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_Interpolator.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_
#define SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
// LNA/src
#include "linalg_typedefs.hpp"
// FEM/INT/src
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
// FEM/MSI/src
#include "cl_MSI_Dof_Type_Enums.hpp"
// GEN/src
#include "GEN_Data_Types.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        class Property;

        class Field_Interpolator
        {
            // tolerance for check
            real mEpsilon = 1e-12;

            // discretization (B-spline) mesh index the field interpolator operates on
            moris_index mDiscretizationMeshIndex = -1;

            // how many fields are to be interpolated
            const uint mNumberOfFields = 0;

            // pointer to space and time interpolation objects
            mtk::Interpolation_Function_Base* mSpaceInterpolation = nullptr;
            mtk::Interpolation_Function_Base* mTimeInterpolation  = nullptr;

            // space and time geometry interpolator
            Geometry_Interpolator* mGeometryInterpolator = nullptr;

            // space, time, and space time number of bases
            uint mNSpaceBases = 0;
            uint mNTimeBases  = 0;
            uint mNFieldBases = 0;

            // space time number of coefficients
            uint mNFieldCoeff = 0;

            // space, time dimensions
            uint mNSpaceDim = 0;
            uint mNTimeDim  = 0;

            // space parametric dimensions
            uint mNSpaceParamDim = 0;

            // parametric point where field is interpolated
            Matrix< DDRMat > mXi;
            Matrix< DDRMat > mTau;

            // matrix of field coefficients uHat
            Matrix< DDRMat > mUHat;

            // dof type
            moris::Cell< MSI::Dof_Type > mDofType;

            // dv type
            moris::Cell< gen::PDV_Type > mDvType;

            // field type
            moris::Cell< mtk::Field_Type > mFieldType;

            // flag for evaluation
            bool mNBuildEval      = true;
            bool mNEval           = true;
            bool mNTransEval      = true;
            bool mdNdxEval        = true;
            bool md2Ndx2Eval      = true;
            bool md3Ndx3Eval      = true;
            bool mdNdtEval        = true;
            bool md2Ndt2Eval      = true;
            bool md2NdxtEval      = true;
            bool mDivOperatorEval = true;

            bool mValEval      = true;
            bool mValTransEval = true;

            bool mGradx1Eval = true;
            bool mGradx2Eval = true;
            bool mGradx3Eval = true;

            bool mGradt1Eval = true;
            bool mGradt2Eval = true;
            bool mGradt3Eval = true;

            bool mGradxtEval = true;

            // storage
            Matrix< DDRMat > mNBuild;
            Matrix< SDRMat > mN;
            Matrix< SDRMat > mNTrans;
            Matrix< DDRMat > mdNdx;
            Matrix< DDRMat > md2Ndx2;
            Matrix< DDRMat > md3Ndx3;
            Matrix< DDRMat > mdNdt;
            Matrix< DDRMat > md2Ndt2;
            Matrix< DDRMat > md2Ndxt;
            Matrix< DDRMat > mDivOperator;

            Matrix< DDRMat > mVal;
            Matrix< DDRMat > mValTrans;

            Matrix< DDRMat > mGradx1;
            Matrix< DDRMat > mGradx2;
            Matrix< DDRMat > mGradx3;

            Matrix< DDRMat > mGradt1;
            Matrix< DDRMat > mGradt2;
            Matrix< DDRMat > mGradt3;

            Matrix< DDRMat > mGradxt;

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
            Field_Interpolator(
                    const uint&                        aNumberOfFields,
                    const mtk::Interpolation_Rule&     aFieldInterpolationRule,
                    Geometry_Interpolator*             aGeometryInterpolator,
                    const moris::Cell< MSI::Dof_Type > aDofType );

            /**
             * constructor
             * @param[ in ] aNumberOfFields           number of interpolated fields
             * @param[ in ] aFieldInterpolationRule   field interpolation rule
             * @param[ in ] aGeometryInterpolator     pointer to geometry interpolator object
             * @param[ in ] aDvType                   dv type for the interpolated fields
             */
            Field_Interpolator(
                    const uint&                    aNumberOfFields,
                    const mtk::Interpolation_Rule& aFieldInterpolationRule,
                    Geometry_Interpolator*         aGeometryInterpolator,
                    const moris::Cell< gen::PDV_Type >  aDvType );

            /**
             * constructor
             * @param[ in ] aNumberOfFields           number of interpolated fields
             * @param[ in ] aFieldInterpolationRule   field interpolation rule
             * @param[ in ] aGeometryInterpolator     pointer to geometry interpolator object
             * @param[ in ] aFieldType                field type for the interpolated fields
             */
            Field_Interpolator(
                    const uint&                          aNumberOfFields,
                    const mtk::Interpolation_Rule&       aFieldInterpolationRule,
                    Geometry_Interpolator*               aGeometryInterpolator,
                    const moris::Cell< mtk::Field_Type > aDvType );

            /**
             * trivial constructor for unit test
             */
            Field_Interpolator(
                    const uint&                        aNumberOfFields,
                    const moris::Cell< MSI::Dof_Type > aDofType )
                    : mNumberOfFields( aNumberOfFields )
                    , mDofType( aDofType )
            {
                mNFieldCoeff = mNumberOfFields;
            }

            /**
             * trivial constructor for unit test
             */
            Field_Interpolator(
                    const uint&                   aNumberOfFields,
                    const moris::Cell< gen::PDV_Type > aDvType )
                    : mNumberOfFields( aNumberOfFields )
                    , mDvType( aDvType )
            {
                mNFieldCoeff = mNumberOfFields;
            }

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
             * reset evaluation flags when coefficients are changed
             */
            void reset_eval_flags_coefficients();

            //------------------------------------------------------------------------------
            /**
             * get dof type
             */
            const moris::Cell< MSI::Dof_Type >&
            get_dof_type() const
            {
                return mDofType;
            }

            //------------------------------------------------------------------------------

            /**
             * @brief Get the discretization mesh index for the DoF type the FI operates on
             *
             * @return const moris_index
             */
            moris_index
            get_discretization_mesh_index() const
            {
                MORIS_ASSERT( mDiscretizationMeshIndex > -1, "Field_Interpolator::get_discretization_mesh_index() - Discretization mesh index not set." );
                return mDiscretizationMeshIndex;
            }

            //------------------------------------------------------------------------------
            /**
             * get dof type
             */
            const moris::Cell< gen::PDV_Type >&
            get_dv_type() const
            {
                return mDvType;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of space dimension
             */
            uint
            get_space_dim() const
            {
                return mNSpaceDim;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of fields
             */
            uint
            get_number_of_fields() const
            {
                return mNumberOfFields;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of space bases
             */
            uint
            get_number_of_space_bases() const
            {
                return mNSpaceBases;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of time bases
             */
            uint
            get_number_of_time_bases() const
            {
                return mNTimeBases;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of space time bases
             */
            uint
            get_number_of_space_time_bases() const
            {
                return mNFieldBases;
            }

            //------------------------------------------------------------------------------
            /**
             * get the number of space time coefficients
             */
            uint
            get_number_of_space_time_coefficients() const
            {
                return mNFieldCoeff;
            }
            //------------------------------------------------------------------------------
            /**
             * get space interpolation order
             */
            mtk::Interpolation_Order
            get_space_interpolation_order()
            {
                return mSpaceInterpolation->get_interpolation_order();
            }

            //------------------------------------------------------------------------------

            void set_discretization_mesh_index( const moris_index aDiscretizationMeshIndex );

            //------------------------------------------------------------------------------
            /**
             * set the parametric point where field is interpolated
             * @param[ in ] aParamPoint evaluation point in space and time
             */
            void set_space_time( const Matrix< DDRMat >& aParamPoint );

            //------------------------------------------------------------------------------
            /**
             * get the parametric point in space where field is interpolated
             */
            const Matrix< DDRMat >&
            get_space() const
            {
                return mXi;
            }

            //------------------------------------------------------------------------------
            /**
             * get the parametric point in time where field is interpolated of tau
             */
            const Matrix< DDRMat >&
            get_time() const
            {
                return mTau;
            }

            //------------------------------------------------------------------------------
            /**
             * set the coefficients of the field uHat
             * @param[ in ] aUHat coefficients
             */
            void set_coeff( const Matrix< DDRMat >& aUHat );

            //------------------------------------------------------------------------------
            /**
             * get the coefficients of the field uHat
             * @param[ out ] mUHat coefficients
             */
            const Matrix< DDRMat >&
            get_coeff() const
            {
                return mUHat;
            }

            //------------------------------------------------------------------------------
            /**
             * return the space time shape functions
             * @param[ out ] shape functions matrix
             *               ( 1 x <number of bases> )
             */
            const Matrix< DDRMat >& NBuild();

            //------------------------------------------------------------------------------
            /**
             * evaluate the space time shape functions
             */
            void eval_NBuild();

            //------------------------------------------------------------------------------
            /**
             * return the N for vector field ( space time shape functions )
             * @param[ out ] ( nNumberOfFields x mNFieldCoeff )
             */
            const Matrix< SDRMat >& N();

            //------------------------------------------------------------------------------
            /**
             * evaluate the N for vector field ( space time shape functions )
             * @param[ out ] ( nNumberOfFields x mNFieldCoeff )
             */
            void eval_N();

            //------------------------------------------------------------------------------
            /**
             * return the transpose of N for vector field ( space time shape functions )
             * @param[ out ] ( nNumberOfFields x mNFieldCoeff )
             */
            const Matrix< SDRMat >& N_trans();

            //------------------------------------------------------------------------------
            /**
             * return the nth order derivatives of the space time shape functions
             * wrt space x
             * @param[ in ]  aDerivativeOrder derivative order
             * @param[ out ] dnNdxn           nth order spatial derivative of the shape functions
             */
            const Matrix< DDRMat >& dnNdxn( const uint& aDerivativeOrder );

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
            const Matrix< DDRMat >& dnNdtn( const uint& aDerivativeOrder );

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
            const Matrix< DDRMat >& d2Ndxt();

            //------------------------------------------------------------------------------
            /**
             * mixed space & time derivative of the space time shape functions
             * wrt time x & t
             */
            void eval_d2Ndxt();

            //------------------------------------------------------------------------------
            /**
             * get the field at given space and time Xi, Tau
             * @param[ out ]          interpolated field
             */
            const Matrix< DDRMat >& val();

            /**
             * evaluates the field at given space and time Xi, Tau
             * @param[ out ]          interpolated field
             */
            void eval_val();

            /**
             * get the transpose of the field at given space and time Xi, Tau
             * @param[ out ]          interpolated field
             */
            const Matrix< DDRMat >& val_trans();

            //------------------------------------------------------------------------------
            /**
             * evaluates the field space derivatives at given space and time evaluation point
             * @param[ in ]  aDerivativeOrder  order of the required derivatives
             * @param[ out ] gradx             space derivatives
             */
            const Matrix< DDRMat >& gradx( const uint& aDerivativeOrder );

            /**
             * evaluates the field space derivatives at given space and time evaluation point
             * @param[ in ]  aDerivativeOrder  order of the required derivatives
             */
            void eval_gradx( const uint& aDerivativeOrder );

            //------------------------------------------------------------------------------
            /**
             * evaluates the field spatial divergence at given space and time evaluation point
             * @return divergence of the field
             */
            moris::real div();

            //------------------------------------------------------------------------------
            /**
             * return the spatial divergence operator at given space and time evaluation point
             * @param[ out ] divergence operator
             */
            const Matrix< DDRMat >& div_operator();

            //------------------------------------------------------------------------------
            /**
             * evaluates the spatial divergence operator at given space and time evaluation point
             */
            void eval_div_operator();

            //------------------------------------------------------------------------------
            /**
             * evaluates the field time derivative at given space and time evaluation point
             * @param[ in ] aDerivativeOrder  order of the required derivative
             * @param[ out ] gradt            time derivatives
             */
            const Matrix< DDRMat >& gradt( const uint& aDerivativeOrder );

            /**
             * evaluates the field time derivative at given space and time evaluation point
             * @param[ in ] aDerivativeOrder  order of the required derivative
             */
            void eval_gradt( const uint& aDerivativeOrder );

            //------------------------------------------------------------------------------
            /**
             * evaluates the mixed field space & time derivative at given space and time evaluation point
             * @param[ out ] gradxt            mixed space time derivatives
             */
            const Matrix< DDRMat >& gradxt();

            /**
             * evaluates the mixed field space & time derivative at given space and time evaluation point
             */
            void eval_gradxt();

            //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_FIELD_INTERPOLATOR_HPP_ */

