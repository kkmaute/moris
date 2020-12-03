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

            uint mCurrentOptAlgInd;    // stores index of current optimization solver

            std::shared_ptr<moris::opt::Problem> mProblem;

            Matrix< DDSMat > mActive; // flag for active/inactive constraints

        public:

            /**
             * Constructor
             */
            Algorithm();

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
            virtual void solve(uint aCurrentOptAlgInd, std::shared_ptr<Problem> aOptProb) = 0;

            /**
             * @brief write restart file with advs as well as upper and lower bounds
             */
            void write_advs_to_file( uint aIterationIndex, const Matrix<DDRMat> aADVs );
        };
    }
}

#endif /* MORIS_CL_OPT_ALGORITHM_HPP_ */
