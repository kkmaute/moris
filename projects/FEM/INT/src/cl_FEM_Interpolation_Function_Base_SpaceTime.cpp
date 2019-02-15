#include "assert.hpp"
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Factory.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function_Base_SpaceTime.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

    Interpolation_Function_Base_SpaceTime::Interpolation_Function_Base_SpaceTime(
            Interpolation_Function_Base* aSpaceInterpolation,
            Interpolation_Function_Base* aTimeInterpolation)
    {
    	mSpaceInterpolation = aSpaceInterpolation;
    	mTimeInterpolation  = aTimeInterpolation;
    }

//------------------------------------------------------------------------------
    Interpolation_Matrix *
	Interpolation_Function_Base_SpaceTime::evalN( Interpolation_Matrix  	& aN,
            									  const Matrix< DDRMat > 	& aXi,
												  const Matrix< DDRMat > 	& aTau)
            {
                // evaluate space shape functions
    			mSpaceInterpolation->eval_N( tSpaceN, aXi );
    	        // evaluate space shape functions
				mTimeInterpolation->eval_N( tTimeN, aTau );

				tN = tSpaceN.matrix_data() * tTimeN.matrix_data();
				aN = tN;
                // return new interpolation function pointer
                return aN;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
