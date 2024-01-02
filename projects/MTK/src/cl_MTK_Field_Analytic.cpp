/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Field_Analytic.cpp
 *
 */

#include "cl_MTK_Field_Analytic.hpp"

#include "fn_dot.hpp"

#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src

namespace moris
{
    namespace mtk
    {
        Field_Analytic::Field_Analytic(
                Vector<Analytic_Field_Function>       aFunction,
                Vector<Analytic_Derivative_Function>  aDerivativeFunction,
                moris::Matrix<DDRMat>              const & aCoefficients,
                Mesh_Pair                                  aMeshPairs,
                uint                               const & aNumberOfFields)
        : Field(aMeshPairs,aNumberOfFields),
          mAnalyticFieldValueFunction(aFunction),
          mAnalyticDerivativeFunction(aDerivativeFunction)
          {
            // check that for each field one function is provided
            MORIS_ERROR( mAnalyticFieldValueFunction.size() == aNumberOfFields,
                    "mtk::Field_Analytic::Field_Analytic - number of value functions needs to equal number of fields.\n");

            MORIS_ERROR( mAnalyticDerivativeFunction.size() == aNumberOfFields,
                    "mtk::Field_Analytic::Field_Analytic - number of derivative functions needs to equal number of fields.\n");

            // check that coefficient vector is column vector
            MORIS_ERROR( aCoefficients.n_cols() == 1,
                    "mtk::Field_Analytic::Field_Analytic - coefficient vector needs to column vector.\n");

            // set coefficient vector
            mCoefficients=aCoefficients;

            // set number of coefficients
            mNumberOfCoefficients=mCoefficients.n_rows();
          }

        // ----------------------------------------------------------------------------------------------

        Field_Analytic::Field_Analytic(
                Analytic_Field_Function             aFunction,
                Analytic_Derivative_Function        aDerivativeFunction,
                moris::Matrix<DDRMat>       const & aCoefficients,
                Mesh_Pair                           aMeshPairs,
                uint                        const & aNumberOfFields)
        : Field(aMeshPairs,aNumberOfFields)
        {
            // assign same function and derivative function to all fields
            for (uint i=0;i<mNumberOfFields;++i)
            {
                mAnalyticFieldValueFunction.push_back ( aFunction );
                mAnalyticDerivativeFunction.push_back ( aDerivativeFunction );
            }

            // check that coefficient vector is column vector
            MORIS_ERROR( aCoefficients.n_cols() == 1,
                    "mtk::Field_Analytic::Field_Analytic - coefficient vector needs to column vector.\n");

            // set coefficient vector
            mCoefficients=aCoefficients;

            // set number of coefficients
            mNumberOfCoefficients=mCoefficients.n_rows();
       }

        // ----------------------------------------------------------------------------------------------

        Field_Analytic::~Field_Analytic()
        {

        }

        // ----------------------------------------------------------------------------------------------

        void Field_Analytic::compute_nodal_values()
        {
            // check that function pointer is set
            for (uint ii=1;ii<mNumberOfFields;++ii)
            {
                MORIS_ASSERT( mAnalyticFieldValueFunction(ii) != nullptr,
                        "mtk::Field_Analytic::compute_nodal_values - value function not set.\n");
            }

            // check that coefficient vector is set
            MORIS_ASSERT( mNumberOfCoefficients > -1,
                    "mtk::Field_Analytic::compute_nodal_values - coefficient vector not set.\n");

            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair.get_interpolation_mesh();

            // make sure that nodal value matrix is properly sized
             mValues.resize( tIPmesh->get_num_nodes(), mNumberOfFields );

            // loop over all nodes
            for (uint tNodeIndex=0;tNodeIndex<mValues.n_rows();++tNodeIndex)
            {
                // loop over all fields
                for (uint tFieldIndex=0;tFieldIndex<mNumberOfFields;++tFieldIndex)
                {
                    mValues( tNodeIndex ,tFieldIndex ) = mAnalyticFieldValueFunction(tFieldIndex)(
                            tIPmesh->get_node_coordinate( tNodeIndex ),
                            mCoefficients);
                }
            }
        }

        // ----------------------------------------------------------------------------------------------

        void Field_Analytic::compute_derivatives_of_field_value(
                Matrix< DDRMat >       & aDerivatives,
                Matrix< IndexMat >     & aCoefIndices,
                uint             const & aNodeIndex,
                uint             const & aFieldIndex)
        {
            // check that function pointer is set
            MORIS_ASSERT( mAnalyticDerivativeFunction(aFieldIndex) != nullptr,
                    "mtk::Field_Analytic::compute_derivatives_of_field_value - derivative funtion not set.\n");

            // check that coefficient vector is set
             MORIS_ASSERT( mNumberOfCoefficients > -1,
                     "mtk::Field_Analytic::compute_nodal_values - coefficient vector not set.\n");

            // assume that analytic function depends on all coefficients
            aCoefIndices.set_size( mNumberOfCoefficients, 1);

            for (uint i=0;i< mCoefficients.n_rows();++i)
            {
                aCoefIndices(i)=i;
            }

            // get interpolation mesh
            mtk::Mesh * tIPmesh = mMeshPair.get_interpolation_mesh();

            // compute derivatives with respect to all coefficients
            mAnalyticDerivativeFunction(aFieldIndex)(
                    tIPmesh->get_node_coordinate( aNodeIndex ),
                    mCoefficients,
                    aDerivatives);

            // check for correct size of derivative vector
            MORIS_ASSERT( (sint)aDerivatives.n_rows() == mNumberOfCoefficients,
                    "mtk::Field_Analytic::compute_derivatives_of_field_value - incorrect dimension of derivative vector.\n");
        }

        // ----------------------------------------------------------------------------------------------
    }
}

