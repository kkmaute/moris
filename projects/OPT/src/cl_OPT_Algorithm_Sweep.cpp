#include "cl_OPT_Algorithm_Sweep.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_sum.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::Algorithm_Sweep(ParameterList aParameterList)
                : Algorithm(aParameterList)
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
            mUpdateObjectives = mParameterList.get<bool>("evaluate_objectives");
            mUpdateConstraints = mParameterList.get<bool>("evaluate_constraints");
            mUpdateObjectiveGradients = mParameterList.get<bool>("evaluate_objective_gradients");
            mUpdateConstraintGradients = mParameterList.get<bool>("evaluate_constraint_gradients");
            mSave = mParameterList.get<bool>("save");
            mPrint = mParameterList.get<bool>("print");
            std::string tFiniteDifferenceType = mParameterList.get<std::string>("finite_difference_type");

            // Set initial compute flags
            mProblem->mUpdateObjectives = this->mUpdateObjectives;
            mProblem->mUpdateConstraints = this->mUpdateConstraints;
            mProblem->mUpdateObjectiveGradients = this->mUpdateObjectiveGradients;
            mProblem->mUpdateConstraintGradients = this->mUpdateConstraintGradients;

            // Get number of ADVs
            uint tNumADVs = mProblem->get_num_advs();

            // Finite differencing
            Matrix<DDRMat> tFiniteDifferenceEpsilons(1, 1);
            string_to_mat(mParameterList.get<std::string>("finite_difference_epsilons"), tFiniteDifferenceEpsilons);
            if (tFiniteDifferenceEpsilons.numel() == 0)
            {
                tFiniteDifferenceEpsilons.set_size(1, 1, 0);
            }
            else
            {
                if (tFiniteDifferenceEpsilons.n_rows() == 1)
                {
                    uint tNumEvals = tFiniteDifferenceEpsilons.n_cols();
                    tFiniteDifferenceEpsilons.resize(tNumADVs, tNumEvals);
                    for (uint tIndex = 1; tIndex < tNumADVs; tIndex++)
                    {
                        tFiniteDifferenceEpsilons({tIndex, tIndex}, {0, tNumEvals - 1}) = tFiniteDifferenceEpsilons({0, 0}, {0, tNumEvals - 1});
                    }
                }
                MORIS_ERROR(tFiniteDifferenceEpsilons.n_rows() == tNumADVs, 
                        "OPT_Algorithm_Sweep: Number of rows in finite_difference_epsilons must match the number of ADVs.");
            }
            uint tTotalEpsilons = tFiniteDifferenceEpsilons.n_cols();

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
                moris::save_matrix_to_hdf5_file(mFileID, "epsilons", tFiniteDifferenceEpsilons, tStatus);
            }

            // Print ADVs/epsilons to be evaluated
            if (mPrint)
            {
                moris::print(tEvaluationPoints, "adv_evaluations");
                moris::print(tFiniteDifferenceEpsilons, "epsilons");
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

                // Get analytical gradients if requested
                if (tFiniteDifferenceType == "none" || tFiniteDifferenceType == "all")
                {
                    mProblem->set_finite_differencing("none");
                    this->evaluate_objective_gradients(tEvaluationName + " analytical");
                    this->evaluate_constraint_gradients(tEvaluationName + " analytical");
                }

                // Loop through epsilons
                for (uint tEpsilonIndex = 0; tEpsilonIndex < tTotalEpsilons; tEpsilonIndex++)
                {

                    // Reset evaluation name with epsilon data
                    tEvaluationName = " eval_" + std::to_string(tEvaluationIndex + 1) + "-" + std::to_string(tTotalEvaluations);
                    if (tFiniteDifferenceType != "none")
                    {
                        tEvaluationName += " epsilon_" + std::to_string(tEpsilonIndex + 1) + "-" + std::to_string(tTotalEpsilons);
                    }

                    // Compute and/or save gradients based on finite differencing requested
                    if (tFiniteDifferenceType == "all")
                    {
                        // Forward
                        mProblem->set_finite_differencing("forward", tFiniteDifferenceEpsilons.get_column(tEpsilonIndex));
                        Matrix<DDRMat> tForwardObjectiveGradient = this->evaluate_objective_gradients(tEvaluationName + " fd_forward");
                        Matrix<DDRMat> tForwardConstraintGradient = this->evaluate_constraint_gradients(tEvaluationName + " fd_forward");

                        // Backward
                        mProblem->set_finite_differencing("backward", tFiniteDifferenceEpsilons.get_column(tEpsilonIndex));
                        Matrix<DDRMat> tBackwardObjectiveGradient = this->evaluate_objective_gradients(tEvaluationName + " fd_backward");
                        Matrix<DDRMat> tBackwardConstraintGradient = this->evaluate_constraint_gradients(tEvaluationName + " fd_backward");

                        // Central
                        this->output_variables((tForwardObjectiveGradient + tBackwardObjectiveGradient) / 2,
                                "objective_gradients" + tEvaluationName + " fd_central");
                        this->output_variables((tForwardConstraintGradient + tBackwardConstraintGradient) / 2,
                                               "constraint_gradients" + tEvaluationName + " fd_central");
                    }
                    else if (tFiniteDifferenceType != "none")
                    {
                        mProblem->set_finite_differencing(tFiniteDifferenceType,
                                                                    tFiniteDifferenceEpsilons.get_column(tEpsilonIndex));
                        this->evaluate_objective_gradients(tEvaluationName);
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
