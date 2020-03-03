// Project header files
#include "cl_OPT_Algorithm_Sweep.hpp" // OPT/src
#include "fn_Parsing_Tools.hpp"
#include "fn_sum.hpp"
// #include "HDF5_Tools.hpp"

// ---------------------------------------------------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Algorithm_Sweep::Algorithm_Sweep() :
                Algorithm()
        {
            // assign default parameter values
            mParameterList.insert("num_evaluations_per_adv", "10");
            mParameterList.insert("custom_adv_evaluations", "");
            mParameterList.insert("include_bounds", true);
            mParameterList.insert("calculate_objectives", true);
            mParameterList.insert("calculate_constraints", true);
            mParameterList.insert("calculate_objective_gradients", true);
            mParameterList.insert("calculate_constraint_gradients", true);
            mParameterList.insert("save", true);
            mParameterList.insert("print", false);
            mParameterList.insert("hdf5_path", "");
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::~Algorithm_Sweep()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_Sweep::solve(Problem* aOptProb )
        {

            // Initialize problem/algorithm
            mProblem = aOptProb; // set the member variable mProblem to aOptProb
            Algorithm::initialize(); // initialize the base class member variables

            // Extract basic sweep parameters
            bool tIncludeBounds = mParameterList.get<bool>("include_bounds");
            bool tUpdateObjectives = mParameterList.get<bool>("calculate_objectives");
            bool tUpdateConstraints = mParameterList.get<bool>("calculate_constraints");
            bool tUpdateObjectiveGradients = mParameterList.get<bool>("calculate_objective_gradients");
            bool tUpdateConstraintGradients = mParameterList.get<bool>("calculate_constraint_gradients");
            bool tSave = mParameterList.get<bool>("save");
            bool tPrint = mParameterList.get<bool>("print");

            // Initialize variables obtained from problem
            uint tNumADVs = mProblem->get_num_advs();
            uint tNumObjectives = mProblem->get_num_objectives();
            uint tNumConstraints = mProblem->get_num_constraints();
            Matrix<DDRMat> tObjectives(tNumObjectives, 1);
            Matrix<DDRMat> tConstraints(tNumConstraints, 1);
            Matrix<DDRMat> tObjectiveGradients(tNumObjectives, tNumADVs);
            Matrix<DDRMat> tConstraintGradients(tNumConstraints, tNumADVs);

            // Get the evaluation points
            Matrix<DDRMat> tEvaluationPoints(1, 1);
            string_to_mat(mParameterList.get<std::string>("custom_adv_evaluations"), tEvaluationPoints);
            uint tTotalEvaluations = tEvaluationPoints.n_cols();
            if (tEvaluationPoints.n_rows() != tNumADVs)
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
                tTotalEvaluations = sum(tNumEvaluations);
                tEvaluationPoints.set_size(tNumADVs, tTotalEvaluations);
                Matrix<DDRMat> tADVs = tLowerBounds;
                Matrix<DDSMat> tCurrentEvaluations(tNumADVs, 1, 1);

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
                    for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                    {
                        if (tCurrentEvaluations(tADVIndex) == tNumEvaluations(tADVIndex))
                        {
                            // Reset this ADV to the lower bound and incremement next ADV
                            tADVs(tADVIndex) = tLowerBounds(tADVIndex);
                            tCurrentEvaluations(tADVIndex) = 1;

                            tADVs(tADVIndex + 1) += (tUpperBounds(tADVIndex) - tLowerBounds(tADVIndex)) / (tNumEvaluations(tADVIndex) + 1 - (2 * tIncludeBounds));
                            tCurrentEvaluations(tADVIndex + 1) += 1;
                        }
                    }
                }
            }

            // Set up file
//            hid_t tFileID = create_hdf5_file(mParameterList.get<std::string>("hdf5_path"));
//            herr_t tStatus = 0;
            
            // Loop through evaluation points
            for (uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++)
            {
                // Set new ADVs
                mProblem->set_advs(tEvaluationPoints.get_column(tEvaluationIndex));

                // Set flags for computation
                mProblem->mUpdateObjectives = tUpdateObjectives;
                mProblem->mUpdateConstraints = tUpdateConstraints;
                mProblem->mUpdateObjectiveGradient = tUpdateObjectiveGradients;
                mProblem->mUpdateConstraintGradient = tUpdateConstraintGradients;

                // Calculate the requested optimization variables
                tObjectives = mProblem->get_objectives();
                tConstraints = mProblem->get_constraints();
                tObjectiveGradients = mProblem->get_objective_gradient();
                tConstraintGradients = mProblem->get_constraint_gradient();

                // Save
                std::string tEvaluationString = "(" + std::to_string(tEvaluationIndex + 1) + "/" + std::to_string(tTotalEvaluations) + ")";
                if (tSave)
                {
//                    save_matrix_to_hdf5_file(tFileID, "ADVs Point #" + tEvaluationString, tEvaluationPoints.get_column(tEvaluationIndex), tStatus);
//                    save_matrix_to_hdf5_file(tFileID, "Objectives Point #" + tEvaluationString, tObjectives, tStatus);
//                    save_matrix_to_hdf5_file(tFileID, "Constraints Point #" + tEvaluationString, tConstraints, tStatus);
//                    save_matrix_to_hdf5_file(tFileID, "Objective Gradients Point #" + tEvaluationString, tObjectiveGradients, tStatus);
//                    save_matrix_to_hdf5_file(tFileID, "Constraint Gradients Point #" + tEvaluationString, tConstraintGradients, tStatus);
                }

                // Print
                if (tPrint)
                {
                    moris::print(tEvaluationPoints.get_column(tEvaluationIndex), tEvaluationString + " ADVs");
                    moris::print(tObjectives, tEvaluationString + " Objectives");
                    moris::print(tConstraints, tEvaluationString + " Constraints");
                    moris::print(tObjectiveGradients, tEvaluationString + " Objective Gradients");
                    moris::print(tConstraintGradients, tEvaluationString + " ConstraintGradients");
                }
            }

            // Close file
//            close_hdf5_file(tFileID);

            aOptProb = mProblem; // update aOptProb
        }
    }
}
