#ifndef MORIS_OPTIMIZATION_CL_OPTPROB_HPP_
#define MORIS_OPTIMIZATION_CL_OPTPROB_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Mat.hpp" // LNA/src

namespace moris
{
    namespace opt
    {
        class OptProb
        {

        private:

            // Declare member variables
            uint mNumMasterAdvs;    // number of advs
            uint mNumCriteria;      // Number of criteria involved
            uint mNumConstraints;   // number of constraints
            uint mNumEqConstraints; // number of equality constraints

            Mat< real > mAbsDesVarVec;    // Abstract Design Variable Vector
            Mat< real > mAbsDesVarVecUp;  // Upper bounds on Abstract Design Variable Vector
            Mat< real > mAbsDesVarVecLow; // Lower bounds on  Abstract Design Variable Vector

            real         mObj;         // objective
            Mat< real > mCon;         // matrix of constraints
            Mat< sint > mTypeCon;     // matrix containing flags for type of constraints
            Mat< real > mOptCriteria; // Set of optimization criteria
            Mat< real > mGradObj;     // derivative of objective w.r.t advs
            Mat< real > mGradCon;     // derivative of constraints w.r.t advs
            Mat< real > mGradCrit;    // derivative of criteria w.r.t advs
            Mat< sint > mActive;      // flag for active/inactive constraints

        public:

            bool mIsOptAnlyt = false; // flag for solving purely analytical problems

            /**
             * Constructor
             */
            OptProb( ) = default;

            /**
             * Constructor
             *
             * @param[in] aNumMasterAdvs The number of master ADVs
             * @param[in] aNumCriteria The number of criteria
             * @param[in] aNumConstraints The number of constraints
             * @param[in] aNumEqConstraints The number of equality constraints
             */
            OptProb(
                    uint aNumMasterAdvs,
                    uint aNumCriteria,
                    uint aNumConstraints,
                    uint aNumEqConstraints );

            /**
             * Destructor
             */
            ~OptProb();

            /**
             * @brief Get number of adv
             */
            uint get_num_adv()
            {
                return mNumMasterAdvs;
            }

            /**
             * @brief Get number of constraints
             */
            uint get_num_con()
            {
                return mNumConstraints;
            }

            /**
             * @brief Get number of equality constraints
             */
            uint get_num_eqcon()
            {
                return mNumEqConstraints;
            }

            /**
             * @brief Get the adv vector
             */
            Mat< real > & get_adv_vec()
            {
                return mAbsDesVarVec;
            }

            /**
             * @brief Get the adv upper bound
             */
            Mat< real > & get_adv_vec_up()
            {
                return mAbsDesVarVecUp;
            }

            /**
             * @brief Get the adv lower bound
             */
            Mat< real > & get_adv_vec_low()
            {
                return mAbsDesVarVecLow;
            }

            /**
             * @brief Get the constraint type vector
             */
            Mat< sint > & get_typecon_vec()
            {
                return mTypeCon;
            }

            /**
             * @brief Define the abstract design variables
             */
            void define_advs();

            /**
             * @brief Define the Optimization criteria
             */
            void define_opt_criteria();

            /**
             * @brief Computes the objective and constraint values
             *
             * @param[in] aIter Optimization iteration number
             * @param[in] aAdvVec The current vector of ADVs
             *
             * @param[out] aObjVal The value of the objective function
             * @param[out] aConVal Vector of value of constraints
             */
            void compute_obj_con(
                    const uint           aIter,
                    const Mat< real > & aAdvVec,
                    real              & aObjVal,
                    Mat< real >       & aConVal );

            /**
             * @brief Compute the criteria based on the forward analysis type
             */
            void compute_criteria();

            /**
             * @brief Computes the gradient of objective and constraint values
             *        w.r.t. abstract design variables
             *
             * @param[in] aActive Vector of active/inactive constraint flags
             *
             * @param[out] aGradObj Derivative of the objective w.r.t. ADVs
             * @param[out] aGradCon Derivative of the constrains w.r.t. ADVs
             */
            void grad_obj_con(
                    const Mat < sint > & aActive,
                    Mat< real >       & aGradObj,
                    Mat< real >       & aGradCon );

            /**
             * @brief Computes the explicit derivative of objective and
             *        constraints w.r.t. ADVs
             *
             * @param[out] aGradObj_dCrit Derivative of the objective w.r.t.
             *             optimization criteria
             * @param[out] aGradCon_dCrit Derivative of the constrains w.r.t.
             *             optimization criteria
             */
            void explicit_dobjcon_dcrit(
                    Mat< real > & aGradObj_dCrit,
                    Mat< real > & aGradCon_dCrit);

            /**
             * @brief Computes the  derivative of the optimization criteria
             *        w.r.t. ADVs
             */
            void dcriteria_ds();

            /**
             * @brief Computes the explicit derivative of objective and
             *        constraints w.r.t. ADVs
             *
             * @param[out] aGradObj_ds Explicit derivative of the objective
             *             w.r.t. ADVs
             * @param[out] aGradCon_ds Derivative of the constrains w.r.t. ADVs
             */
            void explicit_dobjcon_ds(
                    Mat< real > & aGradObj_ds,
                    Mat< real > & aGradCon_ds);

            /**
             * @brief Assembles the overall objective and constraint gradients
             *
             * @param[in] aGradObj_dCrit Derivative of the objective w.r.t.
             *            optimization criteria
             * @param[in] aGradCon_dCrit Derivative of the constrains w.r.t.
             *            optimization criteria
             * @param[in] aGradObj_ds Explicit derivative of the objective
             *            w.r.t. ADVs
             * @param[in] aGradCon_ds Derivative of the constrains w.r.t. ADVs
             */
            void assemble_grad(
                    const Mat< real > & aGradObj_dCrit,
                    const Mat< real > & aGradCon_dCrit,
                    const Mat< real > & aGradObj_ds,
                    const Mat< real > & aGradCon_ds );

            /**
             * @brief Modifies the optimization solution
             */
            void scale_opt_sol();

            /**
             * @brief Update the optimization problem
             */
            void update_opt_prob();
        };
    }
}

#endif /* MORIS_OPTIMIZATION_CL_OPTPROB_HPP_ */
