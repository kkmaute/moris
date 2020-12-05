#include "cl_OPT_Problem.hpp"
#include "op_plus.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Communication_Tools.hpp"
#include "HDF5_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "fn_stringify_matrix.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------
        // Public functions
        // -------------------------------------------------------------------------------------------------------------

        Problem::Problem(
                ParameterList                       & aParameterList,
                std::shared_ptr<Criteria_Interface> & aInterface)
        {
            // Set interface
            mInterface = aInterface;

            // Parameters: restart file name
            mRestartFile = aParameterList.get<std::string>("restart_file");

            // Parameters: finite differencing
            mFiniteDifferenceType = aParameterList.get<std::string>("finite_difference_type");
            string_to_mat(aParameterList.get<std::string>("finite_difference_epsilons"), mFiniteDifferenceEpsilons);
        }

        // -------------------------------------------------------------------------------------------------------------

        Problem::~Problem()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::initialize()
        {
            // Trace optimization problem
            Tracer tTracer( "OPT", "OptProblem", "Initialize" );

            // Initialize ADVs
            mInterface->initialize(mADVs, mLowerBounds, mUpperBounds);

            // Override interface ADVs
            this->override_advs();

            // Log number of optimization variables
            MORIS_LOG_SPEC("Number of optimization variables",mADVs.numel());

            // Check to make sure ADVs are a column vector
            if (mADVs.n_rows() == 1)
            {
                mADVs = trans(mADVs);
            }

            // Set finite difference epsilons knowing number of advs
            this->set_finite_differencing(mFiniteDifferenceType, mFiniteDifferenceEpsilons);

            // Read advs from restart file
            if ( par_rank() == 0 && ! mRestartFile.empty() )
            {
                this->read_restart_file();
            }

            // Initialize constraint types
            mConstraintTypes = this->get_constraint_types();

            MORIS_ERROR(this->check_constraint_order(), "The constraints are not ordered properly (Eq, Ineq)"); //TODO move call to alg

            // Set number of objectives and constraints
            mNumObjectives  = 1;
            mNumConstraints = mConstraintTypes.numel();

            // Log number of objectives and constraints
            MORIS_LOG_SPEC("Number of objectives",mNumObjectives);
            MORIS_LOG_SPEC("Number of objectives",mNumConstraints);
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

        const Matrix<DDRMat>& Problem::get_objectives()
        {
            if (mUpdateObjectives)
            {
                mObjectives = this->compute_objectives();

                MORIS_ASSERT(mObjectives.numel() == mNumObjectives,
                        "Number of objectives is incorrect.\n");

                // log objective
                MORIS_LOG_SPEC( "Objective", ios::stringify_log( mObjectives ) );
            }
            return mObjectives;
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Problem::get_constraints()
        {
            if (mUpdateConstraints)
            {
                mConstraints = this->compute_constraints();

                MORIS_ASSERT(mConstraints.numel() == mNumConstraints,
                        "Number of objectives is incorrect.\n");
               
                // log constraints
                MORIS_LOG_SPEC( "Constraints", ios::stringify_log( mConstraints ) );
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

        void Problem::set_advs( const Matrix<DDRMat> & aNewADVs)
        {
            mADVs = aNewADVs;

            MORIS_ASSERT( mADVs.n_rows() == aNewADVs.n_rows() && mADVs.n_cols() == aNewADVs.n_cols(),
                    "Problem::set_advs - size of ADV vectors does not match.\n");

            mCriteria = mInterface->get_criteria(aNewADVs);

            mInterface->get_dcriteria_dadv();

            // log criteria and ADVs
            MORIS_LOG_SPEC( "Criteria", ios::stringify_log( mCriteria ) );
            if ( mADVs.numel() > 0 )
            {
                MORIS_LOG_SPEC( "MinADV", mADVs.min() );
                MORIS_LOG_SPEC( "MaxADV", mADVs.max() );
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Problem::get_objective_gradients()
        {
            if (mUpdateObjectiveGradients)
            {
                (this->*compute_objective_gradient)();
                mUpdateObjectiveGradients = false;
            }

            return mObjectiveGradient;
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix<DDRMat>& Problem::get_constraint_gradients()
        {
            if (mUpdateConstraintGradients)
            {
                (this->*compute_constraint_gradient)();
                mUpdateConstraintGradients = false;
            }

            return mConstraintGradient;
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::set_finite_differencing(std::string aType, Matrix<DDRMat> aEpsilons)
        {
            // Check input
            if (aEpsilons.numel() == 1)
            {
                aEpsilons.resize(mADVs.length(), 1);
                for (uint tIndex = 1; tIndex < mADVs.length(); tIndex++)
                {
                    aEpsilons(tIndex) = aEpsilons(0);
                }
            }
            MORIS_ERROR(aEpsilons.numel() == mADVs.length(),
                    "OPT_Problem: Number of elements in finite_difference_epsilons must match the number of ADVs.");

            // Assign epsilons
            mFiniteDifferenceEpsilons = aEpsilons;

            // Function pointers
            this->set_finite_differencing(aType);
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::set_finite_differencing(std::string aType)
        {
            // Set gradient function pointers
            switch (aType[0])
            {
                case 'b':
                {
                    mFiniteDifferenceEpsilons = mFiniteDifferenceEpsilons * -1.0;
                }
                case 'f':
                {
                    compute_objective_gradient  = &Problem::compute_objective_gradient_fd_bias;
                    compute_constraint_gradient = &Problem::compute_constraint_gradient_fd_bias;
                    break;
                }
                case 'c':
                {
                    compute_objective_gradient  = &Problem::compute_objective_gradient_fd_central;
                    compute_constraint_gradient = &Problem::compute_constraint_gradient_fd_central;
                    break;
                }
                default:
                {
                    compute_objective_gradient  = &Problem::compute_objective_gradient_analytical;
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
            mObjectiveGradient =
                    this->compute_dobjective_dadv() +
                    this->compute_dobjective_dcriteria() * mInterface->get_dcriteria_dadv();
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_constraint_gradient_analytical()
        {
            mConstraintGradient =
                    this->compute_dconstraint_dadv() +
                    this->compute_dconstraint_dcriteria() * mInterface->get_dcriteria_dadv();
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_objective_gradient_fd_bias()
        {
            // Set perturbed ADVs and objectives
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mObjectiveGradient.set_size(1, mADVs.length());
            real tObjectivePerturbed;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < mADVs.length(); tADVIndex++) 
            {
                // Perturb
                mADVs(tADVIndex) += mFiniteDifferenceEpsilons(tADVIndex);
                mCriteria = mInterface->get_criteria(mADVs);
                tObjectivePerturbed = this->compute_objectives()(0);

                // Biased finite difference
                mObjectiveGradient(0, tADVIndex) =
                        (tObjectivePerturbed - mObjectives(0)) / mFiniteDifferenceEpsilons(tADVIndex);

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_constraint_gradient_fd_bias()
        {
            // Set perturbed ADVs and constraints
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mConstraintGradient.set_size(mConstraints.length(), mADVs.length());
            Matrix<DDRMat> tConstraintsPerturbed;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < mADVs.length(); tADVIndex++)
            {
                // Perturb
                mADVs(tADVIndex) += mFiniteDifferenceEpsilons(tADVIndex);
                mCriteria = mInterface->get_criteria(mADVs);
                tConstraintsPerturbed = this->compute_constraints();

                // Biased finite difference
                for (uint tConstraintIndex = 0; tConstraintIndex < mConstraints.length(); tConstraintIndex++)
                {
                    mConstraintGradient(tConstraintIndex, tADVIndex) =
                            (tConstraintsPerturbed(tConstraintIndex) - mConstraints(tConstraintIndex)) /
                            mFiniteDifferenceEpsilons(tADVIndex);
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
            mObjectiveGradient.set_size(1, mADVs.length());
            real tObjectivePlus;
            real tObjectiveMinus;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < mADVs.length(); tADVIndex++)
            {
                // Perturb forwards
                mADVs(tADVIndex) += mFiniteDifferenceEpsilons(tADVIndex);
                mCriteria = mInterface->get_criteria(mADVs);
                tObjectivePlus = this->compute_objectives()(0);

                // Perturb backwards
                mADVs(tADVIndex) -= 2 * mFiniteDifferenceEpsilons(tADVIndex);
                mCriteria = mInterface->get_criteria(mADVs);
                tObjectiveMinus = this->compute_objectives()(0);

                // Central difference
                mObjectiveGradient(0, tADVIndex) =
                        (tObjectivePlus - tObjectiveMinus) / (2 * mFiniteDifferenceEpsilons(tADVIndex));

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::compute_constraint_gradient_fd_central()
        {
            // Set perturbed ADVs and constraints
            Matrix<DDRMat> tOriginalADVs = mADVs;
            mConstraintGradient.set_size(mConstraints.length(), mADVs.length());
            Matrix<DDRMat> tConstraintsPlus;
            Matrix<DDRMat> tConstraintsMinus;

            // FD each ADV
            for (uint tADVIndex = 0; tADVIndex < mADVs.length(); tADVIndex++)
            {
                // Perturb forwards
                mADVs(tADVIndex) += mFiniteDifferenceEpsilons(tADVIndex);
                mCriteria = mInterface->get_criteria(mADVs);
                tConstraintsPlus = this->compute_constraints();

                // Perturb backwards
                mADVs(tADVIndex) -= 2 * mFiniteDifferenceEpsilons(tADVIndex);
                mCriteria = mInterface->get_criteria(mADVs);
                tConstraintsMinus = this->compute_constraints();

                // Central difference
                for (uint tConstraintIndex = 0; tConstraintIndex < mConstraints.length(); tConstraintIndex++)
                {
                    mConstraintGradient(tConstraintIndex, tADVIndex) =
                            (tConstraintsPlus(tConstraintIndex) - tConstraintsMinus(tConstraintIndex)) / (2.0 * mFiniteDifferenceEpsilons(tADVIndex));
                }

                // Restore ADV
                mADVs(tADVIndex) = tOriginalADVs(tADVIndex);
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void Problem::read_restart_file()
        {
            MORIS_LOG_INFO("Reading ADVs from restart file: %s",mRestartFile.c_str() );

            // Open open restart file
            hid_t tFileID  = open_hdf5_file( mRestartFile );

            // Define matrix in which to read restart ADVs
            Matrix<DDRMat> tRestartADVs;
            Matrix<DDRMat> tRestartUpperBounds;
            Matrix<DDRMat> tRestartLowerBounds;

            // Read ADVS from restart file
            herr_t tStatus = 0;
            load_matrix_from_hdf5_file( tFileID, "ADVs",        tRestartADVs,        tStatus);
            load_matrix_from_hdf5_file( tFileID, "UpperBounds", tRestartUpperBounds, tStatus);
            load_matrix_from_hdf5_file( tFileID, "LowerBounds", tRestartLowerBounds, tStatus);

            // Close restart file
            close_hdf5_file(tFileID);

            // Check for matching sizes
            MORIS_ERROR( tRestartADVs.numel() == mADVs.numel(),
                    "Number of restart ADVS does not match.\n");
            MORIS_ERROR( tRestartUpperBounds.numel() == mUpperBounds.numel(),
                    "Number of restart upper bounds of ADVs does not match.\n");
            MORIS_ERROR( tRestartLowerBounds.numel() == mLowerBounds.numel(),
                    "Number of restart lower bounds of ADVs does not match.\n");

            MORIS_LOG_INFO("Norm of ADV vector - before loading restart %e   after %e.\n",
                    norm(mADVs), norm(tRestartADVs) );

            // Copy restart vectors on member variables
            mADVs        = tRestartADVs;
            mUpperBounds = tRestartUpperBounds;
            mLowerBounds = tRestartLowerBounds;
        }
    }
}
