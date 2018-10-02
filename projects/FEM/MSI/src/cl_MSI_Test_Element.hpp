/*
 * cl_MSI_Test_Element.hpp
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

    class Test_Element : public MSI::Equation_Object
    {
    private:
        Matrix< DDRMat > ( *mFunction )( Matrix< DDRMat > tMyValues, const moris::uint aEquationObjectInd );
    protected:
    public:
        /**
         * constructor
         *
         */
        Test_Element( const moris::Cell< fem::Node_Base * > & aNodeObjs,
                            Matrix< DDRMat > ( *aFunction )(       Matrix< DDRMat > tMyValues,
                                                             const moris::uint   aEquationObjectInd  ) ) : Equation_Object( aNodeObjs )
        {
            mFunction = aFunction;
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Test_Element(){};

//------------------------------------------------------------------------------
        void compute_jacobian_and_residual()
        {
            Matrix< DDRMat > tTMatrix;
            this->build_PADofMap( tTMatrix );

            Matrix< DDSMat > tGlobalIndExtract( 2, 1, 0);
                             tGlobalIndExtract( 1, 0 ) = 1;
            Matrix< DDRMat > tMyValues;
            Matrix< DDRMat > tMyValues2;

            mSolVec->extract_my_values( 2, tGlobalIndExtract, 0, tMyValues );

            tMyValues2 = tTMatrix * tMyValues;

            mResidual = mFunction( tMyValues2, mEqnObjInd );
        };
    };

    }
}

#endif /* SRC_FEM_CL_MSI_TEST_ELEMENT_HPP_ */
