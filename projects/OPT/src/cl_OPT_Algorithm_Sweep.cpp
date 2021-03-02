#include "cl_OPT_Algorithm_Sweep.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_sum.hpp"
#include "HDF5_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::Algorithm_Sweep(ParameterList aParameterList)
                : mIncludeBounds(aParameterList.get<bool>("include_bounds")),
                  mEvaluateObjectives(aParameterList.get<bool>("evaluate_objectives")),
                  mEvaluateConstraints(aParameterList.get<bool>("evaluate_constraints")),
                  mEvaluateObjectiveGradients(aParameterList.get<bool>("evaluate_objective_gradients")),
                  mEvaluateConstraintGradients(aParameterList.get<bool>("evaluate_constraint_gradients")),
                  mSave(aParameterList.get<bool>("save")),
                  mPrint(aParameterList.get<bool>("print")),
                  mFileID(create_hdf5_file(aParameterList.get<std::string>("hdf5_path"))),
                  mFiniteDifferenceType(aParameterList.get<std::string>("finite_difference_type")),
                  mFiniteDifferenceEpsilons(string_to_mat<DDRMat>(aParameterList.get<std::string>("finite_difference_epsilons"))),
                  mNumEvaluations(string_to_mat<DDUMat>(aParameterList.get<std::string>("num_evaluations_per_adv"))),
                  mEvaluationPoints(string_to_mat<DDRMat>(aParameterList.get<std::string>("custom_adv_evaluations"))),
                  mFiniteDifferenceADVs(string_to_mat<DDUMat>(aParameterList.get<std::string>("finite_difference_adv_indices")))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::~Algorithm_Sweep()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_Sweep::solve(
                uint aCurrentOptAlgInd,
                std::shared_ptr<Problem> aOptProb )
        {
            // Trace optimization
            Tracer tTracer( "OptimizationAlgorithm", "Sweep", "Solve" );

            //----------------------------------------------------------------------------------------------------------
            // Sweep not implemented for parallel
            //----------------------------------------------------------------------------------------------------------
            MORIS_ERROR(par_size() == 1, "Sweep algorithm not implemented for parallel.\n");

            //----------------------------------------------------------------------------------------------------------
            // Initialize
            //----------------------------------------------------------------------------------------------------------
            mCurrentOptAlgInd = aCurrentOptAlgInd;  // set index of current optimization algorithm
            mProblem          = aOptProb;           // set the member variable mProblem to aOptProb

            uint tNumADVs = mProblem->get_num_advs();

            // determine with respect to which advs sensitivities are compute by FD; if list not set in input file
            // all sensitivities with respect to all advs are computed
            uint tNumFDadvs = mFiniteDifferenceADVs.numel();
            if ( tNumFDadvs == 0 )
            {
                // Set number of FD-ADVs to number of ADVs
                tNumFDadvs = tNumADVs;

                // fill adv list with all advs
                mFiniteDifferenceADVs.set_size(tNumADVs,1);
                for (uint tIndex=0;tIndex<tNumADVs;++tIndex)
                {
                    mFiniteDifferenceADVs(tIndex)=tIndex;
                }
            }

            // Finite differencing perturbation size
            if (mFiniteDifferenceEpsilons.numel() != 0)
            {
                if (mFiniteDifferenceEpsilons.n_rows() == 1)
                {
                    uint tNumEvals = mFiniteDifferenceEpsilons.n_cols();
                    mFiniteDifferenceEpsilons.resize(tNumADVs, tNumEvals);
                    for (uint tIndex = 1; tIndex < tNumADVs; tIndex++)
                    {
                        mFiniteDifferenceEpsilons({tIndex, tIndex}, {0, tNumEvals - 1}) = mFiniteDifferenceEpsilons({0, 0}, {0, tNumEvals - 1});
                    }
                }
                MORIS_ERROR(mFiniteDifferenceEpsilons.n_rows() == tNumADVs, 
                        "OPT_Algorithm_Sweep: Number of rows in finite_difference_epsilons must match the number of ADVs.");
            }
            uint tTotalEpsilons = mFiniteDifferenceEpsilons.n_cols();

            //----------------------------------------------------------------------------------------------------------
            // Set up evaluation points
            //----------------------------------------------------------------------------------------------------------
            uint tTotalEvaluations = mEvaluationPoints.n_cols();
            if (mEvaluationPoints.numel() == 0)
            {
                // Check user input
                if (mNumEvaluations.numel() == 1) // check for global number of evaluations
                {
                    mNumEvaluations.set_size(tNumADVs, 1, mNumEvaluations(0));
                }
                else // check for number of evaluations given per ADV
                {
                    MORIS_ERROR(mNumEvaluations.numel() == tNumADVs, 
                            "Must give single number of evaluations for all ADVs or one per ADV, or provide custom evaluation points");
                }
                
                // Check lower and upper bounds for equality
                Matrix<DDRMat> tADVs        = mProblem->get_advs();
                Matrix<DDRMat> tLowerBounds = mProblem->get_lower_bounds();
                Matrix<DDRMat> tUpperBounds = mProblem->get_upper_bounds();

                for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                {
                    // If 1 evaluation, set lower and upper bounds to be equal
                    if (mNumEvaluations(tADVIndex) == 1)
                    {
                        tLowerBounds(tADVIndex) = tADVs(tADVIndex);
                        tUpperBounds(tADVIndex) = tADVs(tADVIndex);
                    }

                    // If lower and upper bounds are equal, set 1 evaluation
                    if (tLowerBounds(tADVIndex) == tUpperBounds(tADVIndex))
                    {
                        mNumEvaluations(tADVIndex) = 1;
                    }
                }

                // Change lower bounds for sweep based on parameter and initialize ADvs
                if (!mIncludeBounds)
                {
                    for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                    {
                        tLowerBounds(tADVIndex) += (tUpperBounds(tADVIndex) - tLowerBounds(tADVIndex)) / (mNumEvaluations(tADVIndex) + 1);
                    }
                }
                
                // Set up evaluations
                tTotalEvaluations = 1;
                for (uint ind = 0; ind < mNumEvaluations.numel(); ind++)
                {
                    tTotalEvaluations *= mNumEvaluations(ind);
                }
                mEvaluationPoints.set_size(tNumADVs, tTotalEvaluations);
                tADVs = tLowerBounds;
                Matrix<DDUMat> tCurrentEvaluations(tNumADVs, 1, 0);

                // Construct evaluation points
                for (uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++)
                {
                    if ( tNumADVs > 0 )
                    {
                        // Assign ADVs
                        for (uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++)
                        {
                            mEvaluationPoints(tADVIndex, tEvaluationIndex) = tADVs(tADVIndex);
                        }

                        // Update ADVs
                        tADVs(0) += (tUpperBounds(0) - tLowerBounds(0)) / (mNumEvaluations(0) + 1 - (2 * mIncludeBounds));

                        tCurrentEvaluations(0) += 1;

                        for (uint tADVIndex = 0; tADVIndex < tNumADVs - 1; tADVIndex++)
                        {
                            if (tCurrentEvaluations(tADVIndex) == mNumEvaluations(tADVIndex))
                            {
                                // Reset this ADV to the lower bound and increment next ADV
                                tADVs(tADVIndex) = tLowerBounds(tADVIndex);
                                tCurrentEvaluations(tADVIndex) = 0;

                                tADVs(tADVIndex + 1) += (tUpperBounds(tADVIndex) - tLowerBounds(tADVIndex)) / (mNumEvaluations(tADVIndex) + 1 - (2 * mIncludeBounds));
                                tCurrentEvaluations(tADVIndex + 1) += 1;
                            }
                        }
                    }
                }
            }
            else
            {
                MORIS_ERROR(mEvaluationPoints.n_rows() == tNumADVs, "Number of rows in custom_adv_evaluations must match the number of ADVs.");
            }

            // Open file and write ADVs/epsilons
            herr_t tStatus = 0;
            if (mSave)
            {
                moris::save_matrix_to_hdf5_file(mFileID, "adv_evaluations", mEvaluationPoints, tStatus);
                moris::save_matrix_to_hdf5_file(mFileID, "epsilons", mFiniteDifferenceEpsilons, tStatus);
            }

            // Print ADVs/epsilons to be evaluated
            if (mPrint)
            {
                moris::print(mEvaluationPoints, "adv_evaluations");
                moris::print(mFiniteDifferenceEpsilons, "epsilons");
            }

            //----------------------------------------------------------------------------------------------------------
            // Perform sweep
            //----------------------------------------------------------------------------------------------------------

            // Evaluation string
            std::string tEvaluationName;
            
            // Loop through evaluation points
            for (uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++)
            {
                // Log iteration of optimization
                MORIS_LOG_ITERATION();

                // Set new ADVs
                this->criteria_solve(mEvaluationPoints.get_column(tEvaluationIndex));
                if (mEvaluateObjectiveGradients or mEvaluateConstraintGradients)
                {
                    mProblem->trigger_dcriteria_dadv_solve();
                }

                // Set evaluation name
                tEvaluationName = " eval_" + std::to_string(tEvaluationIndex + 1) + "-" + std::to_string(tTotalEvaluations);

                // Evaluate objectives and constraints
                this->output_objectives_constraints(tEvaluationName);

                // Get analytical gradients if requested
                if (mFiniteDifferenceType == "none" || mFiniteDifferenceType == "all")
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
                    if (mFiniteDifferenceType != "none")
                    {
                        tEvaluationName += " epsilon_" + std::to_string(tEpsilonIndex + 1) + "-" + std::to_string(tTotalEpsilons);
                    }

                    // Compute and/or save gradients based on finite differencing requested
                    if (mFiniteDifferenceType == "all")
                    {
                        // Set indices of ADVs with respect to which sensitivities are computed by FD
                        mProblem->set_finite_differencing_advs( mFiniteDifferenceADVs );

                        // Forward
                        mProblem->set_finite_differencing("forward", mFiniteDifferenceEpsilons.get_column(tEpsilonIndex));
                        Matrix<DDRMat> tForwardObjectiveGradient  = this->evaluate_objective_gradients(tEvaluationName + " fd_forward");
                        Matrix<DDRMat> tForwardConstraintGradient = this->evaluate_constraint_gradients(tEvaluationName + " fd_forward");

                        // Backward
                        mProblem->set_finite_differencing("backward", mFiniteDifferenceEpsilons.get_column(tEpsilonIndex));
                        Matrix<DDRMat> tBackwardObjectiveGradient = this->evaluate_objective_gradients(tEvaluationName + " fd_backward");
                        Matrix<DDRMat> tBackwardConstraintGradient = this->evaluate_constraint_gradients(tEvaluationName + " fd_backward");

                        // Central
                        this->output_variables((tForwardObjectiveGradient + tBackwardObjectiveGradient) / 2,
                                "objective_gradients" + tEvaluationName + " fd_central");
                        this->output_variables((tForwardConstraintGradient + tBackwardConstraintGradient) / 2,
                                               "constraint_gradients" + tEvaluationName + " fd_central");
                    }
                    else if (mFiniteDifferenceType != "none")
                    {
                        // just analytical sensitivity analysis
                        mProblem->set_finite_differencing(mFiniteDifferenceType,
                                                                    mFiniteDifferenceEpsilons.get_column(tEpsilonIndex));
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
            // Output
            if (mEvaluateObjectives)
            {
                Matrix<DDRMat> tObjectives = mProblem->get_objectives();
                this->output_variables(tObjectives, "objectives" + aEvaluationName);
            }
            if (mEvaluateConstraints)
            {
                Matrix<DDRMat> tConstraints = mProblem->get_constraints();
                this->output_variables(tConstraints, "constraints" + aEvaluationName);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Algorithm_Sweep::evaluate_objective_gradients(std::string aEvaluationName)
        {
            if (mEvaluateObjectiveGradients)
            {
                // Calculate the objective gradients
                Matrix<DDRMat> tObjectiveGradients = mProblem->get_objective_gradients();

                // Extract sensitivities with respect to advs specified for FD
                if ( tObjectiveGradients.numel() > mFiniteDifferenceADVs.numel() )
                {
                    // Allocate temporary matrix (needed as adv indices in mFiniteDifferenceADVs may not be sorted
                    Matrix<DDRMat> tObjectiveGradientsExtracted(mFiniteDifferenceADVs.numel(),1);

                    // Extract sensitivities of requested ADVs
                    for (uint tIndex=0;tIndex<mFiniteDifferenceADVs.numel();++tIndex)
                    {
                        tObjectiveGradientsExtracted(tIndex)=tObjectiveGradients(mFiniteDifferenceADVs(tIndex));
                    }

                    // replace full derivative matrix with temporary one
                    tObjectiveGradients=tObjectiveGradientsExtracted;
                }

                // Output
                this->output_variables(tObjectiveGradients, "objective_gradients" + aEvaluationName);

                // Return
                return tObjectiveGradients;
            }
            return {{}};
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Algorithm_Sweep::evaluate_constraint_gradients(std::string aEvaluationName)
        {
            if (mEvaluateConstraintGradients)
            {
                // Calculate the constraint gradients
                Matrix<DDRMat> tConstraintGradients = mProblem->get_constraint_gradients();

                // Extract sensitivities with respect to advs specified for FD
                 if ( tConstraintGradients.n_cols() > mFiniteDifferenceADVs.numel() )
                 {
                     // Number of constraints
                     uint tNumberOfConstraints = tConstraintGradients.n_rows();

                     // Allocate temporary matrix (needed as adv indices in mFiniteDifferenceADVs may not be sorted
                     Matrix<DDRMat> tConstraintGradientsExtracted(tNumberOfConstraints,mFiniteDifferenceADVs.numel(),1);

                     // Extract sensitivities of requested ADVs
                     for (uint tIndex=0;tIndex<mFiniteDifferenceADVs.numel();++tIndex)
                     {
                         for (uint tConst=0;tConst<tNumberOfConstraints;++tConst)
                         {
                             tConstraintGradientsExtracted(tConst,tIndex)=tConstraintGradients(tConst,mFiniteDifferenceADVs(tIndex));
                         }
                     }

                     // replace full derivative matrix with temporary one
                     tConstraintGradients=tConstraintGradientsExtracted;
                 }

                // Output
                this->output_variables(tConstraintGradients, "constraint_gradients" + aEvaluationName);

                // Return
                return tConstraintGradients;
            }
            return {{}};
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_Sweep::output_variables(
                Matrix<DDRMat> aVariables,
                std::string aFullEvaluationName)
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
