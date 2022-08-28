/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * st_GEN_Field_Parameters.hpp
 *
 */

#ifndef MORIS_ST_GEN_FIELD_PARAMETERS_HPP
#define MORIS_ST_GEN_FIELD_PARAMETERS_HPP

#include "cl_Matrix.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * This is a struct used to simplify \ref moris::ge::Field constructors. It contains additional parameters that
         * are used by all fields.
         */
        struct Field_Parameters
        {
            std::string    mName = "";                       //! Name of this field for identification
            Matrix<DDSMat> mNumRefinements = {{}};           //! The number of refinement steps to use for this field
            Matrix<DDSMat> mRefinementMeshIndices = {{}};    //! Indices of meshes to perform refinement on
            sint           mRefinementFunctionIndex = -1;    //! Index of a user-defined refinement function (-1 = default)
            sint           mDiscretizationMeshIndex = -2;    //! Index of a mesh for discretization (-2 = none, -1 = store nodal values)
            real           mDiscretizationLowerBound = -1.0; //! Lower bound for the B-spline coefficients in this field
            real           mDiscretizationUpperBound = 1.0;  //! Upper bound for the B-spline coefficients in this field
        };
    }
}

#endif //MORIS_ST_GEN_FIELD_PARAMETERS_HPP

