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

            uint mCurrentOptAlgInd;                         // stores index of current optimization solver

            std::shared_ptr<moris::opt::Problem> mProblem;  // pointer to problem algorithm operates on

            Matrix< DDSMat > mActive;                       // flag for active/inactive constraints

            bool mRunning = true;

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
              * Sets the new ADVs to the problem and performs a new forward and sensitivity criteria solve in parallel.
              *
              * @param aADVs ADVs, empty if not on proc 0
              */
             void criteria_solve( const Matrix<DDRMat> & aADVs );

             /**
              * Communicates proc 0's running status to other processors so they know when to end.
              */
             void communicate_running_status();

             /**
              * Dummy solve on processors not running optimization algorithm.
              */
             void dummy_solve();

            /**
             * @brief write restart file with advs as well as upper and lower bounds
             */
            void write_advs_to_file( uint aIterationIndex, const Matrix<DDRMat> aADVs );
        };
    }
}

#endif /* MORIS_CL_OPT_ALGORITHM_HPP_ */
