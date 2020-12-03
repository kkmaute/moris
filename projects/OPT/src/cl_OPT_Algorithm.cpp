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

        // -------------------------------------------------------------------------------------------------------------
        
        void Algorithm::write_advs_to_file( uint aIterationIndex, const Matrix<DDRMat> aADVs )
        {
            // Create file name
            std::string tRestartFileName =
                    "ADV_Alg_" +
                    std::to_string(mCurrentOptAlgInd) +
                    "_Iter_" +
                    std::to_string(aIterationIndex) +
                    ".hdf5";

            // Create file
            hid_t tFileID = create_hdf5_file(tRestartFileName);

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
