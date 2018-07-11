    /*
 * cl_tools_Function_Levelset.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: doble
 */

#ifndef SRC_TOOLS_CL_FUNCTION_LEVELSET_HPP_
#define SRC_TOOLS_CL_FUNCTION_LEVELSET_HPP_


#include "cl_Enums.hpp" // TOL/src
#include "cl_Function.hpp" // TOL/src

namespace moris
{
    namespace tools
    {
        class FunctionLevelset : public Function
        {
        public:
            ~FunctionLevelset();

            /*@brief constructs a levelset function, contains a switch for different types of levelset
             * functions
             *
             * @param[in] aFunctionParameters - contains all the necessary paramaters necessary
             * For the supported levelset function types, the structure of aFunctionParameters is detailed below
             * - LEVELSET_SPHERE - (radius, x-center, y-center, z-center)
             *
             */
            FunctionLevelset(moris::Mat<moris::real>          aFunctionParameters,
                    enum FunctionType                aLevelsetType);

        private:
            moris::Mat<moris::real> mFuncParameters;
            FunctionType   mFuncType;



            /* @brief calculates a value for a given coordinate
             *
             * @param[in] aCoords a 1x3 matrix (x,y,z) coordinate
             *
             * @param[out] FuncVal the value of the equation at that point
             */
            moris::real
            compute_func_val(const moris::Mat<moris::real>        & aCoords);

            /* @brief Computes sensitivies dx/dp for each design variable
             *
             * @param[in] aCoords a n_px3 matrix (x,y,z) coordinate (each row corresponds to a design variable dx/dp,dy/dp,dz/dp)
             *
             * @param[out] FuncVal the value of the equation at that point
             */
            moris::Mat<moris::real>
            compute_func_dxdp(const moris::Mat<moris::real>        & aCoords);
        };
    }
}


#endif /* SRC_TOOLS_CL_FUNCTION_LEVELSET_HPP_ */
