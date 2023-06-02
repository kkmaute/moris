/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_T_Matrix_2.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_T_MATRIX_2_HPP_
#define SRC_HMR_CL_HMR_T_MATRIX_2_HPP_

#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_Cell.hpp" //CNT/src
#include "cl_HMR_T_Matrix.hpp" //CNT/src

namespace moris
{
    namespace hmr
    {
        /**
         *  \brief a T-Matrix Generator
         */
//-------------------------------------------------------------------------------

        class T_Matrix_2 : public T_Matrix
        {
            //! ref to Lagrange Mesh
            Lagrange_Mesh_Base * mLagrangeMeshCoarse;

            //! order of Lagrange Mesh
            uint        mLagrangeOrderCoarse;

            //! number of nodes per Lagrange element
            uint        mNumberOfNodesCoarse;

//-------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------

            // constructor
            T_Matrix_2( const Parameters         * aParameters,
                            BSpline_Mesh_Base  * aBSplineMesh,
                            Lagrange_Mesh_Base * aLagrangeMesh,
                            Lagrange_Mesh_Base * aLagrangeMeshCoarse );

//-------------------------------------------------------------------------------

            // destructor
            ~T_Matrix_2();

//-------------------------------------------------------------------------------
            void evaluate( const uint aBSplineMeshIndex,
                           const bool aBool = true);

//-------------------------------------------------------------------------------
       private:
//-------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------
    }
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_T_MATRIX_HPP_ */

