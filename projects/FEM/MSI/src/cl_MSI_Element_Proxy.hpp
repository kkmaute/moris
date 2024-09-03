/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Element_Proxy.hpp
 *
 */

#ifndef SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_
#define SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Node.hpp"    //FEM/INT/src
#include "cl_MSI_Equation_Object.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris::MSI
{
    class Element_Proxy : public MSI::Equation_Object
    {
      private:
        Matrix< DDRMat > ( *mFunction )( Matrix< DDRMat > tMyValues, const moris::uint aEquationObjectInd );

      public:
        /**
         * constructor
         *
         */
        Element_Proxy(
                const Vector< Vector< fem::Node_Base * > > &aNodeObjs,
                Matrix< DDRMat > ( *aFunction )( Matrix< DDRMat > tMyValues,
                        const moris::uint                         aEquationObjectInd ) )
                : Equation_Object( aNodeObjs )
        {
            mFunction = aFunction;
        };

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Element_Proxy() override {};

        //------------------------------------------------------------------------------

        void compute_jacobian() override
        {
        }

        //------------------------------------------------------------------------------

        void compute_residual() override {
            //            Matrix< DDRMat > tTMatrix;
            //            this->build_PADofMap( tTMatrix );
            //
            //            Vector< Matrix< DDRMat > > tMyValues;
            //
            //            mSolVec->extract_my_values( mUniqueAdofList.numel(), mUniqueAdofList, 0, tMyValues );
            //
            //            for( uint Ik = 0; Ik < tMyValues.size(); Ik++ )
            //            {
            //                tMyValues( Ik ) = tTMatrix * tMyValues( Ik );
            //            }
            //
            //            mElementBlock->mResidual.resize( 1 );
            //            mElementBlock->mResidual( 0 ) = mFunction( tMyValues( 0 ), mEqnObjInd );
        };
    };

}    // namespace moris::MSI

#endif /* SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_ */
