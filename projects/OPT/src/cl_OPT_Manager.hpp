#ifndef MORIS_OPTIMIZATION_CL_OPTMANAGER_HPP_
#define MORIS_OPTIMIZATION_CL_OPTMANAGER_HPP_

// C++ header files.
#include <utility>

// MORIS project header files.
#include "cl_Cell.hpp" // CON/src
#include "core.hpp"
#include "cl_OPT_Algorithm_API.hpp" // OPT/src
#include "cl_OPT_Problem.hpp" // OPT/src
#include "cl_Param_List.hpp"

namespace moris
{
    namespace opt
    {
        class Manager
        {
        private:
            Cell< Algorithm_API > mAlgorithms; // Cell of pointers of optimization algorithms
            std::shared_ptr<Problem> mProblem;

        public:
            /**
             * Constructor
             *
             * @param[in] aOptProb Object of type Problem containing relevant
             *            data regarding ADVs, the objective and constraints
             */
            // Manager(Cell< Algorithm_API > aAlgorithms, Problem* aProblem);
            Manager(moris::Cell<moris::Cell<ParameterList>>& tParameterlist);

            /**
             * Destructor
             */
            ~Manager();

            /**
             * Solve the Cell of optimization algorithms
             */
            void perform();

            /**
             * Gets the avs from the problem
             *
             * @return Vector of current advs
             */
            Matrix<DDRMat> get_advs();

            /**
             * Gets the objectives from the problem
             *
             * @return Vector of objectives, aka a single objective value
             */
            Matrix<DDRMat> get_objectives();

        };
    }  // namespace opt
}      // namespace moris

#endif  /* MORIS_OPTIMIZATION_CL_OPTDATA_HPP_ */
