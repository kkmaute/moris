#ifndef MORIS_CL_OPT_MANAGER_HPP_
#define MORIS_CL_OPT_MANAGER_HPP_

#include "cl_Cell.hpp"
#include "cl_Param_List.hpp"
#include "cl_OPT_Algorithm.hpp"
#include "cl_OPT_Problem.hpp"
#include "cl_OPT_Criteria_Interface.hpp"

namespace moris
{
    namespace opt
    {
        class Manager
        {
        private:
            Cell< std::shared_ptr<Algorithm> > mAlgorithms; // Cell of pointers of optimization algorithms
            std::shared_ptr<Problem> mProblem;

        public:

            /**
             * Constructor with both parameter list and a cell of criteria interfaces
             *
             * @param aParameterLists parameter lists for defining an optimization problem
             * @param aCriteriaInterface criteria interfaces, in addition to any specified in the parameter lists
             */
            Manager(const Cell<Cell<ParameterList>>& aParameterLists,
                    Cell<std::shared_ptr<Criteria_Interface>> aInterfaces = Cell<std::shared_ptr<Criteria_Interface>>(0));

            /**
             * Constructor with problem class already created
             *
             * @param aAlgorithmParameterLists parameter lists for defining just the optimization algorithms
             */
            Manager(const Cell<ParameterList>& aAlgorithmParameterLists,
                    std::shared_ptr<Problem> aProblem);

            /**
             * Destructor
             */
            ~Manager();

            /**
             * Solve the Cell of optimization algorithms
             */
            void perform();

            void restart_with_remesh( uint aI );

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
}  // namespace moris

#endif  /* MORIS_CL_OPT_MANAGER_HPP_ */
