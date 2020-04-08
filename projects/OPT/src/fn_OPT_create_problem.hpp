//
// Created by christopherson on 3/4/20.
//

#ifndef MORIS_FN_OPT_CREATE_PROBLEM_HPP
#define MORIS_FN_OPT_CREATE_PROBLEM_HPP

#include "cl_OPT_Problem.hpp"
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        /**
         * Creates an instance of the specified Problem class and returns a shared pointer to it
         *
         * @param aProblemParameterList moris::ParameterList for the optimization problem
         * @param aInterface std::shared_ptr<Interface> pointers to already-built interfaces
         * @return std::shared_ptr<Problem>
         */
        std::shared_ptr<Problem> create_problem(ParameterList aProblemParameterList, std::shared_ptr<Criteria_Interface> aInterface);
    }
}

#endif //MORIS_FN_OPT_CREATE_PROBLEM_HPP
