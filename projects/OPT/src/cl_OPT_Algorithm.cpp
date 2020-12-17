#include "cl_OPT_Algorithm.hpp"
#include "HDF5_Tools.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------

        Algorithm::Algorithm()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        Algorithm::~Algorithm()
        {
        }

        //----------------------------------------------------------------------------------------------------------------------

        void Algorithm::criteria_solve(const Matrix<DDRMat> & aADVs)
        {
            // Log iteration of optimization
            MORIS_LOG_ITERATION();

            // Set ADVs and get criteria
            this->mProblem->set_advs(aADVs);
        }

        //----------------------------------------------------------------------------------------------------------------------

        void Algorithm::communicate_running_status()
        {
            // Sending/receiving status
            Matrix<DDSMat> tSendingStatus = {{mRunning}};
            Matrix<DDSMat> tReceivingStatus(0, 0);

            // Communication list
            Matrix<DDUMat> tCommunicationList(1, 1, 0);
            if (par_rank() == 0)
            {
                // Resize communication list and sending mat
                tCommunicationList.resize(par_size() - 1, 1);
                tSendingStatus.set_size(par_size() - 1, 1, mRunning);

                // Assign communication list
                for (uint tProcessorIndex = 1; tProcessorIndex < (uint)par_size(); tProcessorIndex++)
                {
                    tCommunicationList(tProcessorIndex - 1) = tProcessorIndex;
                }
            }

            // Perform communication
            communicate_scalars(tCommunicationList, tSendingStatus, tReceivingStatus);

            // Assign new status
            if (par_rank() != 0)
            {
                mRunning = tReceivingStatus(0);
            }
        }

        //----------------------------------------------------------------------------------------------------------------------

        void Algorithm::dummy_solve()
        {
            // Communicate that these procs need to start running
            this->communicate_running_status();

            // Create dummy ADVs
            Matrix<DDRMat> tDummyADVs;

            // Keep looping over func/grad calls
            while(mRunning)
            {
                // Call to help out with criteria solve
                this->criteria_solve(tDummyADVs);

                // Communicate running status so these processors know when to exit
                this->communicate_running_status();
            }
        }

        // -------------------------------------------------------------------------------------------------------------
        
        void Algorithm::write_advs_to_file( const Matrix<DDRMat> aADVs )
        {
            // Get iteration from global clock
            uint tOptIter = gLogger.get_opt_iteration();

            // Create file name
            std::string tRestartFileName =
                    "ADV_Alg_" +
                    std::to_string(mCurrentOptAlgInd) +
                    "_Iter_" +
                    std::to_string(tOptIter) +
                    ".hdf5";

            // Create file
            // Note: only processor 0 creates this file; therefore no parallel name extension is used
            hid_t tFileID = create_hdf5_file( tRestartFileName, false );

            // Write advs and upper/lower bounds to file
            herr_t tStatus = 0;

            save_matrix_to_hdf5_file(tFileID, "ADVs", aADVs, tStatus);
            save_matrix_to_hdf5_file(tFileID, "UpperBounds", mProblem->get_upper_bounds(), tStatus);
            save_matrix_to_hdf5_file(tFileID, "LowerBounds", mProblem->get_lower_bounds(), tStatus);

            // Close file
            close_hdf5_file(tFileID);
        }

        // -------------------------------------------------------------------------------------------------------------
    }
}
