#ifndef MORIS_ST_GEN_FIELD_PARAMETERS_HPP
#define MORIS_ST_GEN_FIELD_PARAMETERS_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        struct Field_Parameters
        {
            /**
            * @var mName Name of this field for identification
            * @var mNumRefinements The number of refinement steps to use for this field
            * @var mRefinementMeshIndices Indices of meshes to perform refinement on
            * @var mRefinementFunctionIndex The index of a user-defined refinement function (-1 = {} refinement)
            * @var mBSplineMeshIndex Index of a B-spline mesh for discretization (-2 = none, -1 = store nodal values)
            * @var mBSplineLowerBound The lower bound for the B-spline coefficients describing this field
            * @var mBSplineUpperBound The upper bound for the B-spline coefficients describing this field
             */
            std::string     mName = "";
            Matrix<DDSMat>  mNumRefinements = {{}};
            Matrix<DDSMat>  mRefinementMeshIndices = {{}};
            sint            mRefinementFunctionIndex = -1;
            sint            mBSplineMeshIndex = -2;
            real            mBSplineLowerBound = -1.0;
            real            mBSplineUpperBound = 1.0;
        };
    }
}

#endif //MORIS_ST_GEN_FIELD_PARAMETERS_HPP
