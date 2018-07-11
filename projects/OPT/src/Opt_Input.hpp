#ifndef MORIS_OPTIMIZATION_OPT_INPUT_HPP_
#define MORIS_OPTIMIZATION_OPT_INPUT_HPP_

// MORIS project header files.
#include "cl_Cell.hpp" // CON/src
#include "cl_Opt_Alg_API.hpp" // OPT/src
#include "cl_Opt_Prob.hpp" // OPT/src

namespace moris
{
    namespace opt_input
    {
        // declare global variables such that they can be used within
        // different functions within the input file.
        extern uint        Problem;   // optimization problem type
        extern real        Objective; // objective value based on current advs
        extern Mat< real > ADVs;      // current advs

        /**
         * @brief Call to compute the objective and constraints for unit
         *        test purposes
         *
         * @param[in]  aProblem Optimization test problem Id.
         * @param[out] aObjective Value of the objective function
         * @param[out] aAdvVec Vector of value of the current advs
         */
        void create_solve_opt_problem(
                uint          aProblem,
                real        & aObjective,
                Mat< real > & aAdvVec );

        /**
         * @brief Get the solution strategy for the Optimization Problem
         *
         * @return Vector of optimization algorithms
         */
        Cell< opt::OptAlgAPI > define_opt_sol_strategy( );

        /**
         * @brief Define the abstract design variables
         *
         * @param[out] aAbsDesVarVec Vector of the ADV
         * @param[out] aAbsDesVarVecUp Upper bound of ADV
         * @param[out] aAbsDesVarVecLow Lower bound of ADV
         */
        void define_advs(
                Mat< real > & aAbsDesVarVec,
                Mat< real > & aAbsDesVarVecUp,
                Mat< real > & aAbsDesVarVecLow );

        /**
         * @brief Define the Optimization criteria
         *
         * @param[out] aOptCriteria Vector of the optimization criteria
         */
        void define_opt_criteria( Mat< real > & aOptCriteria );

        /**
         * @brief get the analytical criteria
         *
         * @param[in]  aAbsDesVarVec Vector of the abstract design variables
         *
         * @param[out] aCriteria Vector of the set of criteria
         */
        void get_criteria(
                const Mat< real > & aAbsDesVarVec,
                Mat< real >       & aCriteria );

        /**
         * @brief get the derivative of analytical criteria
         *
         * @param[in]  aAbsDesVarVec Vector of the abstract design variables
         *
         * @param[out] aGradCriteria Vector of the gradient of criteria
         */
        void get_dcriteria_ds(
                const Mat< real > & aAbsDesVarVec,
                Mat< real >       & aGradCriteria );

        /**
         * @brief get the objective and constraint values
         *
         * @param[in]  aIter Optimization iteration
         * @param[in]  aAbsDesVarVec Vector of the abstract design variables
         * @param[in]  aCriteria Vector of the set of criteria
         *
         * @param[out] aObjective Value of the objective function
         * @param[out] aConstraints Vector of value of the constraints
         * @param[out] aTypeCon Flag for type of constraints
         */
        void get_obj_con(
                const uint          aIter,
                const Mat< real > & aAbsDesVarVec,
                const Mat< real > & aCriteria,
                real              & aObjective,
                Mat< real >       & aConstraints,
                Mat< sint >       & aTypeCon );

        /**
         * @brief get the derivaitve objective and constraint values w.r.t. adv
         *
         * @param[in]  aAbsDesVarVec Vector of the abstract design variables
         * @param[in]  aCriteria     Vector of the set of criteria
         *
         * @param[out] aDObjective_Ds   Derivative of objective w.r.t. adv
         * @param[out] aDConstraints_Ds Derivative of constraints w.r.t. adv
         */
        void get_dobjcon_ds(
                const Mat< real > & aAbsDesVarVec,
                const Mat< real > & aCriteria,
                Mat< real >       & aDObjective_Ds,
                Mat< real >       & aDConstraints_Ds );

        /**
         * @brief get the derivaitve objective and constraint values w.r.t. criteria
         *
         * @param[in]  aAbsDesVarVec Vector of the abstract design variables
         * @param[in]  aCriteria Vector of the set of criteria
         *
         * @param[out] aDObjective_DCrit Derivative of objective w.r.t. criteria
         * @param[out] aDConstraints_DCrit Derivative of constraints w.r.t. criteria
         */
        void get_dobjcon_dcrit(
                Mat< real >   aAbsDesVarVec,
                Mat< real >   aCriteria,
                Mat< real > & aDObjective_DCrit,
                Mat< real > & aDConstraints_DCrit );
    }
}

#endif /* MORIS_OPTIMIZATION_OPT_INPUT_HPP_ */
