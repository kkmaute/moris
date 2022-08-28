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

#include "typedefs.hpp"           //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Node.hpp"         //FEM/INT/src
#include "cl_MSI_Equation_Object.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace MSI
    {
        class Element_Proxy : public MSI::Equation_Object
        {
        private:
            Matrix< DDRMat > ( *mFunction )( Matrix< DDRMat > tMyValues, const moris::uint aEquationObjectInd );

            fem::Set * mElementBlock;

        protected:

        public:
            /**
             * constructor
             *
             */
            Element_Proxy( const moris::Cell< moris::Cell< fem::Node_Base * > > & aNodeObjs,
                           Matrix< DDRMat > ( *aFunction )(       Matrix< DDRMat > tMyValues,
                                                                  const moris::uint      aEquationObjectInd  ) ) : Equation_Object( aNodeObjs )
            {
                mFunction = aFunction;
            };

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~Element_Proxy(){};

            //------------------------------------------------------------------------------

            void compute_jacobian()
            {

            }

            //------------------------------------------------------------------------------

            void compute_residual()
            {
                //            Matrix< DDRMat > tTMatrix;
                //            this->build_PADofMap( tTMatrix );
                //
                //            moris::Cell< Matrix< DDRMat > > tMyValues;
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

    }
}

#endif /* SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_ */

