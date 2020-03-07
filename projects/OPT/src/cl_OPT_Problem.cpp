// Project header files
#include "cl_OPT_Problem.hpp" // OPT/src
#include "op_plus.hpp"
#include "fn_norm.hpp"
#include "fn_OPT_create_interface.hpp"

extern moris::Logger gLogger;

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------
        // Public functions
        // -------------------------------------------------------------------------------------------------------------

        Problem::Problem(ParameterList aParameterList)
        {
            // Create interface
            mInterface = create_interface(aParameterList);

            // Parameters: finite differencing
            set_objective_finite_differencing(aParameterList.get<std::string>("objective_finite_difference_type"),
                                               aParameterList.get<real>("objective_finite_difference_epsilon"));
            set_constraint_finite_differencing(aParameterList.get<std::string>("constraint_finite_difference_type"),
                                               aParameterList.get<real>("constraint_finite_difference_epsilon"));

        }

        // -------------------------------------------------------------------------------------------------------------

        Problem::~Problem()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::initialize()
        {
            // Initialize ADVs
            mADVs = mInterface->initialize_advs(); // get some ADVs from the interface
            this->override_advs(); // user can override the interface ADVs
            mInterface->begin_new_analysis(mADVs); // potentially new ADVs set and passed back to interface to compute criteria

            // Get the criteria at the first step
            mCriteria = mInterface->get_criteria();

            // Initialize bounds and constraints
            mLowerBounds = mInterface->get_lower_adv_bounds();
            mUpperBounds = mInterface->get_upper_adv_bounds();
            mConstraintTypes = this->get_constraint_types();
            MORIS_ERROR(this->check_constraint_order(), "The constraints are not ordered properly (Eq, Ineq)"); //TODO move call to alg
        }

        // -------------------------------------------------------------------------------------------------------------

        bool Problem::check_constraint_order()
        {
            // Checks that the constraint order is  1. equality constraints   (typ=0)
            //                                      2. inequality constraints (typ=1)

            int tSwitches = 1;
            for (uint tConstraintIndex = 0; tConstraintIndex < mConstraintTypes.numel(); tConstraintIndex++)
            {
                if (mConstraintTypes(tConstraintIndex) == tSwitches)
                {
                    tSwitches -= 1;
                }
            }
            return (tSwitches >= 0);
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem::get_objectives()
        {
            if (mUpdateObjectives)
            {
                mObjectives = this->compute_objectives();
            }
            return mObjectives;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem::get_constraints()
        {
            if (mUpdateConstraints)
            {
                mConstraints = this->compute_constraints();
            }
            return mConstraints;
        }

        // -------------------------------------------------------------------------------------------------------------

        uint Problem::get_num_equality_constraints()
        {
            uint tNumEqualityConstraints = 0;

            // Loop over the constraint types matrix
            for (uint tConstraintIndex = 0; tConstraintIndex < mConstraintTypes.numel(); tConstraintIndex++)
            {
                if (mConstraintTypes(tConstraintIndex) == 0)
                {
                    tNumEqualityConstraints++;
                }
            }

            return tNumEqualityConstraints;
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::set_advs(Matrix<DDRMat> aNewADVs)
        {
            if (norm(aNewADVs - mADVs) < mADVNormTolerance)
            {
                return;
            }

            mADVs = aNewADVs;
            mInterface->begin_new_analysis(aNewADVs);
            mCriteria = mInterface->get_criteria();
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem::get_objective_gradients()
        {
            if (mUpdateObjectiveGradients)
            {
                (this->*compute_objective_gradient)();
                mUpdateObjectiveGradients = false;
            }

            return mObjectiveGradient;
        }

        // -------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Problem::get_constraint_gradients()
        {
            if (mUpdateConstraintGradients)
            {
                (this->*compute_constraint_gradient)();
                mUpdateConstraintGradients = false;
            }

            return mConstraintGradient;
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::set_objective_finite_differencing(std::string aType, real aEpsilon)
        {
            // Set new epsilon if it isn't zero
            if (aEpsilon != 0.0)
            {
                mObjectiveFiniteDifferenceEpsilon = aEpsilon;
            }
            mObjectiveFiniteDifferenceEpsilon = std::abs(mObjectiveFiniteDifferenceEpsilon);

            // Set objective gradient function pointer
            switch (aType[0])
            {
                case('b'):
                {
                    mObjectiveFiniteDifferenceEpsilon *= -1;
                }
                case('f'):
                {
                    compute_objective_gradient = &Problem::compute_objective_gradient_fd_bias;
                    break;
                }
                case('c'):
                {
                    compute_objective_gradient = &Problem::compute_objective_gradient_fd_central;
                    break;
                }
                default:
                {
                    compute_objective_gradient = &Problem::compute_objective_gradient_analytical;
                }
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::set_constraint_finite_differencing(std::string aType, real aEpsilon)
        {
            // Set new epsilon if it isn't zero
            if (aEpsilon != 0.0)
            {
                mConstraintFiniteDifferenceEpsilon = aEpsilon;
            }
            mConstraintFiniteDifferenceEpsilon = std::abs(mConstraintFiniteDifferenceEpsilon);

            // Set constraint gradient function pointer
            switch (aType[0])
            {
                case('b'):
                {
                    mConstraintFiniteDifferenceEpsilon *= -1;
                }
                case('f'):
                {
                    compute_constraint_gradient = &Problem::compute_constraint_gradient_fd_bias;
                    break;
                }
                case('c'):
                {
                    compute_constraint_gradient = &Problem::compute_constraint_gradient_fd_central;
                    break;
                }
                default:
                {
                    compute_constraint_gradient = &Problem::compute_constraint_gradient_analytical;
                }
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::scale_solution()
        {
            // TODO Need to decide on a framework for the scaling of the solution
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::update_problem()
        {
            // TODO Need to decide on a framework for the update of the problem
        }

        // -------------------------------------------------------------------------------------------------------------
        // Private: possible functions for computing gradients (analytical, forward/backward/central finite difference)
        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_objective_gradient_analytical()
        {
            mObjectiveGradient = this->compute_dobjective_dadv()
                                 + this->compute_dobjective_dcriteria() * mInterface->get_dcriteria_dadv();
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_constraint_gradient_analytical()
        {
            mConstraintGradient = this->compute_dconstraint_dadv()
                                  + this->compute_dconstraint_dcriteria() * mInterface->get_dcriteria_dadv();
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_objective_gradient_fd_bias()
        {
            // Set perturbed ADVs and objectives
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mObjectiveGradient.set_size(1, this->get_num_advs());
            real tObjectivePerturbed;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < this->get_num_advs(); tADVIndex++) 
            {
                // Perturb
                mADVs(tADVIndex) += mObjectiveFiniteDifferenceEpsilon;
                mInterface->begin_new_analysis(mADVs);
                mCriteria = mInterface->get_criteria();
                tObjectivePerturbed = this->compute_objectives()(0);

                // Biased finite difference
                mObjectiveGradient(0, tADVIndex) =
                        (tObjectivePerturbed - mObjectives(0)) / mObjectiveFiniteDifferenceEpsilon;

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_constraint_gradient_fd_bias()
        {
            // Set perturbed ADVs and constraints
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mConstraintGradient.set_size(this->get_num_constraints(), this->get_num_advs());
            Matrix<DDRMat> tConstraintsPerturbed;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < this->get_num_advs(); tADVIndex++)
            {
                // Perturb
                mADVs(tADVIndex) += mConstraintFiniteDifferenceEpsilon;
                mInterface->begin_new_analysis(mADVs);
                mCriteria = mInterface->get_criteria();
                tConstraintsPerturbed = this->compute_constraints();

                // Biased finite difference
                for (uint tConstraintIndex = 0; tConstraintIndex < this->get_num_constraints(); tConstraintIndex++)
                {
                    mConstraintGradient(tConstraintIndex, tADVIndex)
                            = (tConstraintsPerturbed(tConstraintIndex) - mConstraints(tConstraintIndex)) / mConstraintFiniteDifferenceEpsilon;
                }

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_objective_gradient_fd_central()
        {
            // Set perturbed ADVs and objectives
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mObjectiveGradient.set_size(1, this->get_num_advs());
            real tObjectivePlus;
            real tObjectiveMinus;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < this->get_num_advs(); tADVIndex++)
            {
                // Perturb forwards
                mADVs(tADVIndex) += mObjectiveFiniteDifferenceEpsilon;
                mInterface->begin_new_analysis(mADVs);
                mCriteria = mInterface->get_criteria();
                tObjectivePlus = this->compute_objectives()(0);

                // Perturb backwards
                mADVs(tADVIndex) -= 2 * mObjectiveFiniteDifferenceEpsilon;
                mInterface->begin_new_analysis(mADVs);
                tObjectiveMinus = this->compute_objectives()(0);

                // Central difference
                mObjectiveGradient(0, tADVIndex) =
                        (tObjectivePlus - tObjectiveMinus) / (2 * mObjectiveFiniteDifferenceEpsilon);

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_constraint_gradient_fd_central()
        {
            // Set perturbed ADVs and constraints
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mConstraintGradient.set_size(this->get_num_constraints(), this->get_num_advs());
            Matrix<DDRMat> tConstraintsPlus;
            Matrix<DDRMat> tConstraintsMinus;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < this->get_num_advs(); tADVIndex++)
            {
                // Perturb forwards
                mADVs(tADVIndex) += mConstraintFiniteDifferenceEpsilon;
                mInterface->begin_new_analysis(mADVs);
                mCriteria = mInterface->get_criteria();
                tConstraintsPlus = this->compute_constraints();

                // Perturb backwards
                mADVs(tADVIndex) -= 2 * mConstraintFiniteDifferenceEpsilon;
                mInterface->begin_new_analysis(mADVs);
                tConstraintsMinus = this->compute_constraints();

                // Central difference
                for (uint tConstraintIndex = 0; tConstraintIndex < this->get_num_constraints(); tConstraintIndex++)
                {
                    mConstraintGradient(tConstraintIndex, tADVIndex)
                    = (tConstraintsPlus(tConstraintIndex) - tConstraintsMinus(tConstraintIndex)) / (2 * mConstraintFiniteDifferenceEpsilon);
                }

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

    }
}