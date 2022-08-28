/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Git_info.hpp
 *
 */

#ifndef PROJECTS_MRS_IOS_SRC_CL_GIT_INFO_HPP_
#define PROJECTS_MRS_IOS_SRC_CL_GIT_INFO_HPP_

#include <string>

namespace moris
{
    // -----------------------------------------------------------------------------

    class git_info
    {
        private:
            std::string mMORIS_GIT_BRANCH;
            std::string mMORIS_GIT_HASH;

        public:
            git_info();

            std::string get_git_branch();

            std::string get_git_hash();
    };
}

#endif

