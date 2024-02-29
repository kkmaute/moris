/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_User_Defined_Refinement.cpp
 *
 */

#include <string>
#include <iostream>

#include "cl_HMR_Element.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
             * It operates on the Lagrange mesh. If an element is flagged,
             * HMR remembers the next active ancestor on the B-Spline mesh
             * and flags it as well.
             *
             * @param[ in ] aElement                 : Lagrange element the flagging is performed on
             * @param[ in ] aElementLocalNodeValues  : Cell of local node Values for the field of interest
             * @param[ in ] aParameters              : relevant parameters
             *
             * The function returns the following settings
             *
             * -1 : do not flag this element for refinement
             *  0 : keep this element by flagging this element
             *  1 : flag this element for refinement
             */
            int
            user_defined_refinement(
                          Element                  * aElement,
                    const Cell< Matrix< DDRMat > > & aElementLocalValues,
                          ParameterList            & aParameters )
            {
                /**
                 * Basic Element Functions

                       // return the moris index of this element
                       aElement->get_index();

                       // return the moris id if this element
                       aElement->get_id();

                       // return the proc owner of this element
                       aElement->get_owner();
                 */

                /**
                 * functions concerning minumum refinement criterion
                 *
                       // returns the minumum level
                       aElement->get_min_refinement_level

                       // forces the minumum level to a value
                       aElement->set_min_refinement_level( uint aMinumumLevel )

                       // only takes the value if the stored one is lower
                       aElement->update_min_refinement_level( uint aMinumumLevel )
                 */

                /**
                 * Grab values from parameter list
                 */
                // max level
                uint tMaxLevel   = aParameters.get< sint >("max_level");

                // minimal value
                real tLowerBound = aParameters.get< real >("lower_bound");

                // maximal
                real tUpperBound = aParameters.get< real >("upper_bound");

                // example of use:
                if (    aElementLocalValues( 0 ).max() >= tLowerBound
                     && aElementLocalValues( 0 ).min() <= tUpperBound )
                {
                    if( aElement->get_level() < tMaxLevel )
                    {
                            return  1;  // refine
                    }
                    else
                    {
                        return 0 ; // keep
                    }
                }
                else
                {
                    return -1;  // do not refine
                }
            }

//------------------------------------------------------------------------------
        }
    }
//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

