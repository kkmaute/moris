/*
 * cl_MSI_Element_Proxy.hpp
 *
 *  Created on: Oct 01, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_
#define SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_

#include "typedefs.hpp"           //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Node.hpp"         //FEM/INT/src
#include "cl_MSI_Equation_Object.hpp"
#include "cl_Vector.hpp"

namespace moris
{
    namespace MSI
    {

    class Element_Proxy : public MSI::Equation_Object
    {
    private:
        Matrix< DDRMat > ( *mFunction )( Matrix< DDRMat > tMyValues, const moris::uint aEquationObjectInd );
    protected:
    public:
        /**
         * constructor
         *
         */
        Element_Proxy( const moris::Cell< fem::Node_Base * > & aNodeObjs,
                            Matrix< DDRMat > ( *aFunction )(       Matrix< DDRMat > tMyValues,
                                                             const moris::uint   aEquationObjectInd  ) ) : Equation_Object( aNodeObjs )
        {
            mFunction = aFunction;
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Element_Proxy(){};

//------------------------------------------------------------------------------

        void compute_jacobian_and_residual()
        {
            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap( tTMatrix );

            Matrix< DDRMat > tMyValues;

            mSolVec->extract_my_values( mUniqueAdofList.length(), mUniqueAdofList, 0, tMyValues );

            tMyValues = tTMatrix * tMyValues;

            mResidual = mFunction( tMyValues, mEqnObjInd );
        };
    };

    }
}

#endif /* SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_ */
