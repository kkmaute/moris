/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Debug.hpp
 *
 */

#ifndef SRC_TOOLS_CL_DEBUG_HPP_
#define SRC_TOOLS_CL_DEBUG_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
namespace moris
{
    class Debug
    {
    public:

        /**
         * Get duplicates of a coordinate list
         *
         * @param[in]  aCoord         .... Coordinate list with x,y,z
         * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
         *
         */
       static Matrix < DDUMat >
       duplicate_row_check(Matrix< DDRMat >  & aCoord);

       /**
        * Get duplicates of an Id list
        *
        * @param[in]  aId            .... Id list of coordinates
        * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
        *
        */
       static Matrix < DDUMat >
           duplicate_row_check(Matrix < DDUMat >  & aId);

       /**
        * Get duplicates of two Id lists
        *
        * @param[in]  aId1           .... Id list of coordinates
        * @param[in]  aId2           .... Another Id list of coordinates
        * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
        *
        */
       static Matrix < DDUMat >
           duplicate_row_check(Matrix < DDUMat >  & aId1,
                               Matrix < DDUMat >  & aId2);

       /**
        * Generates a duplicate list, but returns only entries that have problems between two Id lists
        *
        * @param[in]  aId1         .... Id list of coordinates
        * @param[in]  aId2         .... Another Id list of coordinates
        * @param[out] problem_list .... Shows the problem list [Position(i) Position(j)]
        *
        */
       static Matrix < DDUMat >
       duplicate_row_check_problems(Matrix < DDUMat >  & aId1,
                                        Matrix < DDUMat >  & aId2);

        static Matrix < DDUMat >
        duplicate_col_check(Matrix < DDUMat >  & aId);

    };

}

#endif /* SRC_TOOLS_CL_DEBUG_HPP_ */

