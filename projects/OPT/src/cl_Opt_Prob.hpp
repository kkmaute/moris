#ifndef MORIS_OPTIMIZATION_CL_OPTPROB_HPP_
#define MORIS_OPTIMIZATION_CL_OPTPROB_HPP_

// MORIS project header files.
#include "core.hpp"
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"

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

            Matrix< DDRMat >  mAbsDesVarVec;    // Abstract Design Variable Vector
            Matrix< DDRMat >  mAbsDesVarVecUp;  // Upper bounds on Abstract Design Variable Vector
            Matrix< DDRMat >  mAbsDesVarVecLow; // Lower bounds on  Abstract Design Variable Vector

            real         mObj;         // objective
            Matrix< DDRMat >  mCon;         // matrix of constraints
            Matrix< DDSMat >  mTypeCon;     // matrix containing flags for type of constraints
            Matrix< DDRMat >  mOptCriteria; // Set of optimization criteria
            Matrix< DDRMat >  mGradObj;     // derivative of objective w.r.t advs
            Matrix< DDRMat >  mGradCon;     // derivative of constraints w.r.t advs
            Matrix< DDRMat >  mGradCrit;    // derivative of criteria w.r.t advs
            Matrix< DDSMat >  mActive;      // flag for active/inactive constraints

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
            Matrix< DDRMat >  & get_adv_vec()
            {
                return mAbsDesVarVec;
            }

            /**
             * @brief Get the adv upper bound
             */
            Matrix< DDRMat >  & get_adv_vec_up()
            {
                return mAbsDesVarVecUp;
            }

            /**
             * @brief Get the adv lower bound
             */
            Matrix< DDRMat >  & get_adv_vec_low()
            {
                return mAbsDesVarVecLow;
            }

            /**
             * @brief Get the constraint type vector
             */
            Matrix< DDSMat >  & get_typecon_vec()
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
                    const Matrix< DDRMat >  & aAdvVec,
                    real              & aObjVal,
                    Matrix< DDRMat >        & aConVal );

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
                    const Matrix< DDSMat >  & aActive,
                    Matrix< DDRMat >        & aGradObj,
                    Matrix< DDRMat >        & aGradCon );

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
                    Matrix< DDRMat >  & aGradObj_dCrit,
                    Matrix< DDRMat >  & aGradCon_dCrit);

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
                    Matrix< DDRMat >  & aGradObj_ds,
                    Matrix< DDRMat >  & aGradCon_ds);

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
                    const Matrix< DDRMat >  & aGradObj_dCrit,
                    const Matrix< DDRMat >  & aGradCon_dCrit,
                    const Matrix< DDRMat >  & aGradObj_ds,
                    const Matrix< DDRMat >  & aGradCon_ds );

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
