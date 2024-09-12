/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_OPT_create_algorithm.hpp
 *
 */

#ifndef MORIS_FN_OPT_CREATE_ALGORITHM_HPP
#define MORIS_FN_OPT_CREATE_ALGORITHM_HPP

#include "cl_OPT_Algorithm.hpp"
#include "cl_Parameter_List.hpp"

namespace moris::opt
{
    /**
     * Creates an instance of the specified Algorithm class and returns a shared pointer to it.
     *
     * @param aAlgorithmParameterList Parameter list for an OPT algorithm
     * @return Algorithm class
     */
    std::shared_ptr< Algorithm > create_algorithm( const Parameter_List& aAlgorithmParameterList );
    }

#endif //MORIS_FN_OPT_CREATE_ALGORITHM_HPP

