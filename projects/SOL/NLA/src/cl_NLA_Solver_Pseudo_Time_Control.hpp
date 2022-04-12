/*
 * cl_NLA_Solver_Load_Control.hpp
 */

#ifndef SRC_FEM_CL_NLA_SOLVER_PSEUDO_TIME_CONTROL_HPP_
#define SRC_FEM_CL_NLA_SOLVER_PSEUDO_TIME_CONTROL_HPP_

#include "cl_Param_List.hpp"

#include "cl_SOL_Enums.hpp"

namespace moris
{
    class Solver_Interface;

    namespace sol
    {
        class SOL_Warehouse;
        class Dist_Vector;
    }    // namespace sol

    namespace NLA
    {
        class Solver_Pseudo_Time_Control
        {
          private:
            /// relaxation strategy
            sol::SolverPseudoTimeControlType mTimeStepStrategy;

            /// time step counter
            uint mTimeStepCounter = 1;

            /// Constant time step size
            real mConstantStepSize = 1.0;

            /// Initial time step size used in very first time step
            real mInitialStepSize = 1e18;

            /// pre-factor for time step index-based increase of time step
            real mTimeStepIndexFactor = 0.0;

            /// exponent for time step index-based increase
            real mTimeStepExponent = 1.0;

            /// pre-factor for residual-based increase of time step
            real mResidualFactor = 0.0;

            /// exponent for residual-based increase
            real mResidualExponent = 1.0;

            /// Comsol parameter (laminar: 20, 2D turbulent: 25, 3D turbulent: 30)
            real mComsolParameter_1 = 20.0;

            /// Comsol parameter (laminar: 40, 2D turbulent: 50, 3D turbulent: 60)
            real mComsolParameter_2 = 40.0;

            /// required pseudo time step size needed for convergence
            real mRequiredStepSize = 1000.0;

            /// maximum time step
            real mMaxStepSize = 1e18;

            /// maximum number of time steps
            uint mMaxNumTimeSteps = 1;

            /// required relative residual drop to update "previous" solution and time step size
            real mRelativeResidualDropThreshold;

            /// time step size used once maximum number of step has been reached
            /// for negative values, time step size defined by strategy is used
            real mSteadyStateStepSize = -1.0;

            /// total pseudo time
            real mTotalTime = 0.0;

            /// time offset for outputting pseudo time solutions
            real mTimeOffSet = 0.0;

            /// pseudo time for outputting pseudo time solutions
            real mOutputTime = 0.0;

            // solver interface
            Solver_Interface* mSolverInterface = nullptr;

            // vector for "previous solution"
            sol::Dist_Vector* mPreviousSolution = nullptr;

            // vector with true previous solution (set by time solver)
            sol::Dist_Vector* mOldPrevTimeStepVector = nullptr;

            // true time frames (set by time solver)
            Matrix< DDRMat > mTimeFrameCurrent;
            Matrix< DDRMat > mTimeFramePrevious;

          public:
            Solver_Pseudo_Time_Control(
                    ParameterList&      aParameterListNonlinearSolver,
                    sol::Dist_Vector*   aCurrentSolution,
                    Solver_Interface*   aSolverInterface,
                    sol::SOL_Warehouse* aSolverWarehouse );

            ~Solver_Pseudo_Time_Control();

            /*
             * get initial step size and set convergence status
             */
            bool get_initial_step_size( real& aInitialStepSize );

            /*
             * compute the new time step size and check for convergence
             */
            bool compute_time_step_size(
                    const real&       aRefNorm,
                    const real&       aResNorm,
                    sol::Dist_Vector* aCurrentSolution,
                    real&             aTimeStep,
                    real&             aTotalTime );
        };
    }    // namespace NLA
}    // namespace moris

#endif /* SRC_FEM_CL_NLA_SOLVER_PSEUDO_TIME_CONTROL_HPP_ */
