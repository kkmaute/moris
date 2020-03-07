// Project header files
#include "cl_OPT_Algorithm_Sweep.hpp" // OPT/src
#include "fn_Parsing_Tools.hpp"
#include "fn_sum.hpp"
#include "HDF5_Tools.hpp"

// ---------------------------------------------------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Algorithm_Sweep::Algorithm_Sweep() : Algorithm()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::~Algorithm_Sweep()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_Sweep::solve(std::shared_ptr<Problem> aOptProb )
        {
            //----------------------------------------------------------------------------------------------------------
            // Initialize
            //----------------------------------------------------------------------------------------------------------
            mProblem = aOptProb; // set the member variable mProblem to aOptProb

            // Extract basic sweep parameters
            bool tIncludeBounds = mParameterList.get<bool>("include_bounds");
            mUpdateObjectives = mParameterList.get<bool>("compute_objectives");
            mUpdateConstraints = mParameterList.get<bool>("compute_constraints");
            mUpdateObjectiveGradients = mParameterList.get<bool>("calculate_objective_gradients");
            mUpdateConstraintGradients = mParameterList.get<bool>("calculate_constraint_gradients");
            mSave = mParameterList.get<bool>("save");
            mPrint = mParameterList.get<bool>("print");
            std::string tObjectiveFiniteDifferenceType = mParameterList.get<std::string>("objective_finite_difference_type");
            std::string tConstraintFiniteDifferenceType = mParameterList.get<std::string>("constraint_finite_difference_type");
            
            // Finite difference epsilons
            Matrix<DDRMat> tObjectiveFiniteDifferenceEpsilons(1, 1);
            Matrix<DDRMat> tConstraintFiniteDifferenceEpsilons(1, 1);
            string_to_mat(mParameterList.get<std::string>("objective_finite_difference_epsilons"), tObjectiveFiniteDifferenceEpsilons);
            string_to_mat(mParameterList.get<std::string>("constraint_finite_difference_epsilons"), tConstraintFiniteDifferenceEpsilons);
            if (tObjectiveFiniteDifferenceEpsilons.numel() == 0)
            {
                tObjectiveFiniteDifferenceEpsilons.set_size(1, 1, 0);
            }
            if (tConstraintFiniteDifferenceEpsilons.numel() == 0)
            {
                tConstraintFiniteDifferenceEpsilons.set_size(1, 1, 0);
            }
            MORIS_ERROR(tObjectiveFiniteDifferenceEpsilons.numel() == tConstraintFiniteDifferenceEpsilons.numel(), 
                    "objective_finite_difference_epsilons and constraint_finite_difference_epsilons must have the same number of elements");
            uint tTotalEpsilons = tObjectiveFiniteDifferenceEpsilons.numel();

            // Set initial compute flags
            mProblem->mUpdateObjectives = this->mUpdateObjectives;
            mProblem->mUpdateConstraints = this->mUpdateConstraints;
            mProblem->mUpdateObjectiveGradients = this->mUpdateObjectiveGradients;
            mProblem->mUpdateConstraintGradients = this->mUpdateConstraintGradients;

            // Get number of ADVs
            uint tNumADVs = mProblem->get_num_advs();

            //----------------------------------------------------------------------------------------------------------
            // Set up evaluation points
            //----------------------------------------------------------------------------------------------------------
            Matrix<DDRMat> tEvaluationPoints(1, 1);
            string_to_mat(mParameterList.get<std::string>("custom_adv_evaluations"), tEvaluationPoints);
            uint tTotalEvaluations = tEvaluationPoints.n_cols();
            if (tEvaluationPoints.numel() == 0)
            {
                // Set up based on number of evaluations per adv
                Matrix<DDSMat> tNumEvaluations(1, 1);
                std::string tStringEvaluations = mParameterList.get<std::string>("num_evaluations_per_adv");
                string_to_mat(tStringEvaluations, tNumEvaluations);
                
                // Check user input
                if (tNumEvaluations.numel() == 1) // check for global number of evaluations
                {
                    tNumEvaluations.set_size(tNumADVs, 1, tNumEvaluations(0));
                }
                else // check for number of evaluations given per ADV
                {
                    MORIS_ERROR(tNumEvaluations.numel() == tNumADVs, "Must give single number of evaluations for all ADVs or one per ADV, or provide custom evaluation points");
                }
                
                // Check lower and upper bounds for equality
                Matrix<DDRMat> tLowerBounds = mProblem->get_lower_bounds();
                Matrix<DDRMat> tUpperBounds = mProblem->get_upper_bounds();
                for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                {
                    if (tLowerBounds(tADVIndex) == tUpperBounds(tADVIndex))
                    {
                        tNumEvaluations(tADVIndex) = 1;
                    }
                }

                // Change lower bounds for sweep based on parameter and initialize ADvs
                if (!tIncludeBounds)
                {
                    for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                    {
                        tLowerBounds(tADVIndex) += (tUpperBounds(tADVIndex) - tLowerBounds(tADVIndex)) / (tNumEvaluations(tADVIndex) + 1);
                    }
                }
                
                // Set up evaluations
                tTotalEvaluations = 1;
                for (uint ind = 0; ind < tNumEvaluations.numel(); ind++)
                {
                    tTotalEvaluations *= tNumEvaluations(ind);
                }
                tEvaluationPoints.set_size(tNumADVs, tTotalEvaluations);
                Matrix<DDRMat> tADVs = tLowerBounds;
                Matrix<DDSMat> tCurrentEvaluations(tNumADVs, 1, 0);

                // Construct evaluation points
                for (uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++)
                {
                    // Assign ADVs
                    for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                    {
                        tEvaluationPoints(tADVIndex, tEvaluationIndex) = tADVs(tADVIndex);
                    }

                    // Update ADVs
                    tADVs(0) += (tUpperBounds(0) - tLowerBounds(0)) / (tNumEvaluations(0) + 1 - (2 * tIncludeBounds));
                    tCurrentEvaluations(0) += 1;
                    for (uint tADVIndex = 0; tADVIndex < tNumADVs - 1; tADVIndex++)
                    {
                        if (tCurrentEvaluations(tADVIndex) == tNumEvaluations(tADVIndex))
                        {
                            // Reset this ADV to the lower bound and incremement next ADV
                            tADVs(tADVIndex) = tLowerBounds(tADVIndex);
                            tCurrentEvaluations(tADVIndex) = 0;

                            tADVs(tADVIndex + 1) += (tUpperBounds(tADVIndex) - tLowerBounds(tADVIndex)) / (tNumEvaluations(tADVIndex) + 1 - (2 * tIncludeBounds));
                            tCurrentEvaluations(tADVIndex + 1) += 1;
                        }
                    }
                }
            }
            else
            {
                MORIS_ERROR(tEvaluationPoints.n_rows() == tNumADVs, "Number of rows in custom_adv_evaluations must match the number of ADVs.");
            }


            // Open file and write ADVs/epsilons
            mFileID = create_hdf5_file(mParameterList.get<std::string>("hdf5_path"));
            herr_t tStatus = 0;
            if (mSave)
            {
                moris::save_matrix_to_hdf5_file(mFileID, "adv_evaluations", tEvaluationPoints, tStatus);
                moris::save_matrix_to_hdf5_file(mFileID, "objective_epsilons", tObjectiveFiniteDifferenceEpsilons, tStatus);
                moris::save_matrix_to_hdf5_file(mFileID, "constraint_epsilons", tConstraintFiniteDifferenceEpsilons, tStatus);
            }

            // Print ADVs/epsilons to be evaluated
            if (mPrint)
            {
                moris::print(tEvaluationPoints, "adv_evaluations");
                moris::print(tObjectiveFiniteDifferenceEpsilons, "objective_epsilons");
                moris::print(tConstraintFiniteDifferenceEpsilons, "constraint_epsilons");
            }

            //----------------------------------------------------------------------------------------------------------
            // Perform sweep
            //----------------------------------------------------------------------------------------------------------

            // Evaluation string
            std::string tEvaluationName;
            
            // Loop through evaluation points
            for (uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++)
            {
                // Set new ADVs
                mProblem->set_advs(tEvaluationPoints.get_column(tEvaluationIndex));

                // Set evaluation name
                tEvaluationName = " eval_" + std::to_string(tEvaluationIndex + 1) + "-" + std::to_string(tTotalEvaluations);

                // Evaluate objectives and constraints
                this->output_objectives_constraints(tEvaluationName);

                // Loop through epsilons
                for (uint tEpsilonIndex = 0; tEpsilonIndex < tTotalEpsilons; tEpsilonIndex++)
                {
                    // Add epsilon data
                    if ((tObjectiveFiniteDifferenceType != "none") || (tConstraintFiniteDifferenceType != "none"))
                    {
                        tEvaluationName += " epsilon_" + std::to_string(tEpsilonIndex + 1) + "-" + std::to_string(tTotalEpsilons);
                    }

                    // Compute and save objective gradients based on finite differencing requested
                    if (tObjectiveFiniteDifferenceType == "all")
                    {
                        // Analytical
                        mProblem->set_objective_finite_differencing("none",
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        this->evaluate_objective_gradients(tEvaluationName + " analytical");

                        // Forward
                        mProblem->set_objective_finite_differencing("forward",
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        Matrix<DDRMat> tForwardGradient = this->evaluate_objective_gradients(
                                tEvaluationName + " fd_forward");

                        // Backward
                        mProblem->set_objective_finite_differencing("backward",
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        Matrix<DDRMat> tBackwardGradient = this->evaluate_objective_gradients(
                                tEvaluationName + " fd_backward");

                        // Central
                        this->output_variables((tForwardGradient + tBackwardGradient) / 2,
                                "objective_gradients" + tEvaluationName + " fd_central");
                    }
                    else
                    {
                        mProblem->set_objective_finite_differencing(tObjectiveFiniteDifferenceType,
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        this->evaluate_objective_gradients(tEvaluationName);
                    }

                    // Compute and save constraint gradients based on finite differencing requested
                    if (tConstraintFiniteDifferenceType == "all")
                    {
                        // Analytical
                        mProblem->set_constraint_finite_differencing("none",
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        this->evaluate_constraint_gradients(tEvaluationName + " analytical");

                        // Forward
                        mProblem->set_constraint_finite_differencing("forward",
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        Matrix<DDRMat> tForwardGradient = this->evaluate_constraint_gradients(
                                tEvaluationName + " fd_forward");

                        // Backward
                        mProblem->set_constraint_finite_differencing("backward",
                                                                    tObjectiveFiniteDifferenceEpsilons(tEpsilonIndex));
                        Matrix<DDRMat> tBackwardGradient = this->evaluate_constraint_gradients(
                                tEvaluationName + " fd_backward");

                        // Central
                        this->output_variables((tForwardGradient + tBackwardGradient) / 2,
                                               "constraint_gradients" + tEvaluationName + " fd_central");
                    }
                    else
                    {
                        mProblem->set_constraint_finite_differencing(tConstraintFiniteDifferenceType,
                                                                     tConstraintFiniteDifferenceEpsilons(tEpsilonIndex));
                        this->evaluate_constraint_gradients(tEvaluationName);
                    }
                }
            }

            // Close file
            close_hdf5_file(mFileID);

            // Update problem
            aOptProb = mProblem;
        }
        
        //--------------------------------------------------------------------------------------------------------------
        
        void Algorithm_Sweep::output_objectives_constraints(std::string aEvaluationName)
        {
            // Set flags for computation
            mProblem->mUpdateObjectives = this->mUpdateObjectives;
            mProblem->mUpdateConstraints = this->mUpdateConstraints;

            // Calculate the objectives and constraints
            Matrix<DDRMat> tObjectives = mProblem->get_objectives();
            Matrix<DDRMat> tConstraints = mProblem->get_constraints();

            // Output
            this->output_variables(tObjectives, "objectives" + aEvaluationName);
            this->output_variables(tConstraints, "constraints" + aEvaluationName);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Algorithm_Sweep::evaluate_objective_gradients(std::string aEvaluationName)
        {
            // Set flags for computation
            mProblem->mUpdateObjectiveGradients = this->mUpdateObjectiveGradients;

            // Calculate the objective gradients
            Matrix<DDRMat> tObjectiveGradients = mProblem->get_objective_gradients();

            // Output
            this->output_variables(tObjectiveGradients, "objective_gradients" + aEvaluationName);

            // Return
            return tObjectiveGradients;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Algorithm_Sweep::evaluate_constraint_gradients(std::string aEvaluationName)
        {
            // Set flags for computation
            mProblem->mUpdateConstraintGradients = this->mUpdateConstraintGradients;

            // Calculate the constraint gradients
            Matrix<DDRMat> tConstraintGradients = mProblem->get_constraint_gradients();

            // Output
            this->output_variables(tConstraintGradients, "constraint_gradients" + aEvaluationName);

            // Return
            return tConstraintGradients;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_Sweep::output_variables(Matrix<DDRMat> aVariables, std::string aFullEvaluationName)
        {
            // Write status
            herr_t tStatus = 0;

            // Save
            if (mSave)
            {
                moris::save_matrix_to_hdf5_file(mFileID, aFullEvaluationName, aVariables, tStatus);
            }

            // Print
            if (mPrint)
            {
                moris::print(aVariables, aFullEvaluationName);
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}
