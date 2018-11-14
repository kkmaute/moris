#include <string>
#include <iostream>
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_HMR_Element.hpp"

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
    namespace moris
    {
        namespace hmr
        {
//------------------------------------------------------------------------------

            /**
             * This function returns true if an element fits the criterion.
             *
             * @param[ in ] aElementLocalNodeValues  : Node Values for the field of interest
             * @param[ in ] aParameters              : Cell of relevant parameters
             */
            bool
            user_defined_refinement(
                    const Element                  * aElement,
                    const Cell< Matrix< DDRMat > > & aElementLocalValues,
                    ParameterList                  & aParameters )
            {
                // max level
                uint tMaxLevel   = aParameters.get< sint >("max_level");

                // minimal value
                real tLowerBound =  aParameters.get< real >("lower_bound");

                return  aElementLocalValues( 0 ).max() >= tLowerBound && aElement->get_level() < tMaxLevel;
            }

//------------------------------------------------------------------------------
        }
    }
//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif
