#include <string>

#include "cl_HMR_State.hpp"

namespace moris
{
    namespace hmr
    {
//--------------------------------------------------------------------------------
        class Arguments
        {
            std::string mParameterPath = "";
            std::string mHdf5InputPath  = "";
            std::string mHdf5OutputPath = "";
            std::string mExodusPath    = "";
            bool        mTensorFlag    = false;
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
            get_hdf5_input_path() const
            {
                return mHdf5InputPath;
            }

//---------------------------------------------------------------------------------

            /**
             * return the output path
             */
            const std::string &
            get_hdf5_output_path() const
            {
                return mHdf5OutputPath;
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
        };
//---------------------------------------------------------------------------------
    }
}
