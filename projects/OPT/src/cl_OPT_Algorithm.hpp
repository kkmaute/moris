#ifndef MORIS_CL_OPT_ALGORITHM_HPP_
#define MORIS_CL_OPT_ALGORITHM_HPP_

#include "cl_Param_List.hpp"
#include "cl_OPT_Problem.hpp"

namespace moris
{
    namespace opt
    {
        class Algorithm
        {
        protected:
            std::shared_ptr<moris::opt::Problem> mProblem;
            Matrix< DDSMat > mActive; // flag for active/inactive constraints
            ParameterList mParameterList;
            bool mPrint = false;

        public:

            /**
             * Constructor
             */
            Algorithm(ParameterList aParameterList);

            /**
             * Destructor
             */
            virtual ~Algorithm();

            /**
             * @brief Calls the derived optimization algorithm
             *
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            virtual void solve(std::shared_ptr<Problem> aOptProb) = 0;
        };
    }
}

#endif /* MORIS_CL_OPT_ALGORITHM_HPP_ */
