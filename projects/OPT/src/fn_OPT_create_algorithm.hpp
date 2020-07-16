#ifndef MORIS_FN_OPT_CREATE_ALGORITHM_HPP
#define MORIS_FN_OPT_CREATE_ALGORITHM_HPP

#include "cl_OPT_Algorithm.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        /**
         * Creates an instance of the specified Algorithm class and returns a shared pointer to it.
         *
         * @param aAlgorithmParameterList Parameter list for an OPT algorithm
         * @return Algorithm class
         */
        std::shared_ptr<Algorithm> create_algorithm(ParameterList aAlgorithmParameterList);
    }
}

#endif //MORIS_FN_OPT_CREATE_ALGORITHM_HPP
