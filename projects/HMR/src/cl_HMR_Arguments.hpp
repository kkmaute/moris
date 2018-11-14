#include <string>

#include "cl_HMR_State.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------
        class Arguments
        {
            std::string mParameterPath  = "";
            std::string mDatabaseInputPath  = "";
            std::string mDatabaseOutputPath = "";
            std::string mExodusPath     = "";
            std::string mLastStepPath   = "";
            std::string mBinaryPath     = "";
            std::string mCoeffsPath     = "";
            State       mState;
            double      mTimestep = 0.0;
//--------------------------------------------------------------------------------
        public:
//--------------------------------------------------------------------------------

            Arguments(
                    int  & argc,
                    char * argv[] );

//---------------------------------------------------------------------------------

            void
            print_usage();

//---------------------------------------------------------------------------------

            void
            print_help();

//---------------------------------------------------------------------------------

            /**
             * return the run state of the executable
             */
            State
            get_state() const
            {
                return mState;
            }

//---------------------------------------------------------------------------------

            /**
             * return the parameter path
             */
            const std::string &
            get_parameter_path() const
            {
                return mParameterPath;
            }

//---------------------------------------------------------------------------------

            /**
             * return the input path
             */
            const std::string &
            get_database_input_path() const
            {
                return mDatabaseInputPath;
            }

//---------------------------------------------------------------------------------

            /**
             * return the path of the last step
             */
            const std::string &
            get_last_step_path() const
            {
                return mLastStepPath;
            }

//---------------------------------------------------------------------------------

            /**
             * return the output path
             */
            const std::string &
            get_database_output_path() const
            {
                return mDatabaseOutputPath;
            }

 //---------------------------------------------------------------------------------


            /**
             * return the exodus output path
             */
            const std::string &
            get_exodus_output_path() const
            {
                return mExodusPath;
            }

//---------------------------------------------------------------------------------

            /**
             * return the timestep variable
             */
            double
            get_timestep() const
            {
                return mTimestep;
            }

//---------------------------------------------------------------------------------

            const std::string &
            get_binary_path() const
            {
                return mBinaryPath;
            }

 //---------------------------------------------------------------------------------

            const std::string &
            get_coeffs_path() const
            {
                return mCoeffsPath;
            }

//---------------------------------------------------------------------------------
        };
//---------------------------------------------------------------------------------
    }
}
