/*
 * cl_Interpolation.hpp
 *
 *  Created on: Mar 20, 2017
 *      Author: doble
 */

#ifndef SRC_TOOLS_CL_INTERPOLATION_HPP_
#define SRC_TOOLS_CL_INTERPOLATION_HPP_

#include "cl_Mat.hpp" // LNA/src

namespace moris
{
    class Interpolation
    {
    public:

        //http://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch16.d/IFEM.Ch16.pdf
        /*
         * Linear interpolation of a value given a location
         *
         * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
         * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
         */

       static moris::Mat<moris::real>
        linear_interpolation_location(const moris::Mat<moris::real>  & aInterpVars,
                                      const moris::Mat<moris::real>  & aLocation);

       /*
        * Linear interpolation to a location given a specific value
        */
       static moris::Mat<moris::real>
        linear_interpolation_value(const moris::Mat<moris::real>  & aInterpVars,
                                   const moris::Mat<moris::real>  & aValue);

       /*
        * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
        * Requires 1 local coordinate
        *
        * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
        * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
        */

       static moris::Mat<moris::real>
        bilinear_interpolation(const moris::Mat<moris::real>  & aInterpVars,
                               const moris::Mat<moris::real>  & aValue);


       /*
        * Linear interpolation based on a local coordinate (aLclCoords) based on interpolation vars (aInterpVars)
        * Requires 1 local coordinate
        *
        * @param[in] aInterpVars - Interpolation Vars (x,y,z are treated as independent interpolation variables
        * @param[in] aLclCoords  - Local coordinates to interpolate to (a point at the center of edge has {{0}}
        */
       static moris::Mat<moris::real>
        trilinear_interpolation(const moris::Mat<moris::real>  & aInterpVars,
                                const moris::Mat<moris::real>  & aValue);
    };

}



#endif /* SRC_TOOLS_CL_INTERPOLATION_HPP_ */
