#ifndef MORIS_OPTIMIZATION_CL_OPTALG_HPP_
#define MORIS_OPTIMIZATION_CL_OPTALG_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Param_List.hpp" // CON/src
#include "cl_OPT_Problem.hpp" // OPT/src

namespace moris
{
    namespace opt
    {
        class Algorithm
        {
        public:
            std::shared_ptr<moris::opt::Problem> mProblem; // Object of type optimization problem
            Matrix< DDSMat > mActive; // flag for active/inactive constraints
            ParameterList mParameterList; // The Algorithm specific parameter list

        public:
            bool mPrint = false;

            /**
             * Constructor
             */
            Algorithm();

            /**
             * @brief virtual copy constructor through cloning
             */
            virtual Algorithm*
            clone() const = 0;

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

            /**
             * @brief Initialize the member variables
             */
            void initialize();

            /**
             * @brief Accessor to set a value in the parameter list of Algorithm
             *
             * @param[in] aKey Key corresponding to the mapped value that
             *            needs to be accessed
             */
            ParameterListTypes& set_param( char const* aKey)
            {
                return mParameterList(aKey);
            }
        };
    }  // namespace opt
}      // namespace moris

#endif /* MORIS_OPTIMIZATION_CL_OPTALG_HPP_ */
