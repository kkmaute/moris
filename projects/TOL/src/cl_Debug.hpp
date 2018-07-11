/*
 * cl_Debug.hpp
 *
 *  Created on: Apr 10, 2017
 *      Author: gleim
 */

#ifndef SRC_TOOLS_CL_DEBUG_HPP_
#define SRC_TOOLS_CL_DEBUG_HPP_

#include "cl_Mat.hpp" // LNA/src
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
       static moris::Mat<moris::uint>
           duplicate_row_check(moris::Mat<moris::real>  & aCoord);

       /**
        * Get duplicates of an Id list
        *
        * @param[in]  aId            .... Id list of coordinates
        * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
        *
        */
       static moris::Mat<moris::uint>
           duplicate_row_check(moris::Mat<moris::uint>  & aId);

       /**
        * Get duplicates of two Id lists
        *
        * @param[in]  aId1           .... Id list of coordinates
        * @param[in]  aId2           .... Another Id list of coordinates
        * @param[out] duplicate_list .... Shows the duplicates [Position(i) Position(j)]
        *
        */
       static moris::Mat<moris::uint>
           duplicate_row_check(moris::Mat<moris::uint>  & aId1,
                                                   moris::Mat<moris::uint>  & aId2);

       /**
        * Generates a duplicate list, but returns only entries that have problems between two Id lists
        *
        * @param[in]  aId1         .... Id list of coordinates
        * @param[in]  aId2         .... Another Id list of coordinates
        * @param[out] problem_list .... Shows the problem list [Position(i) Position(j)]
        *
        */
       static moris::Mat<moris::uint>
           duplicate_row_check_problems(moris::Mat<moris::uint>  & aId1,
                                                            moris::Mat<moris::uint>  & aId2);

        static moris::Mat<moris::uint>
        duplicate_col_check(moris::Mat<moris::uint>  & aId);

    };

}



#endif /* SRC_TOOLS_CL_DEBUG_HPP_ */
