// Project header files
#include "cl_Opt_Alg_Sweep.hpp" // OPT/src

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        OptAlgSweep::OptAlgSweep() :
            OptAlg()
        {
            // assign default parameter values
            mParameterList.insert( "max_its", 10 );
        }

        //----------------------------------------------------------------------

        OptAlgSweep::~OptAlgSweep()
        {
        }

        //----------------------------------------------------------------------

        void OptAlgSweep::solve( OptProb & aOptProb )
        {

            mOptProb = aOptProb; // set the member variable mOptProb to aOptProb

            OptAlg::initialize(); // initialize the base class member variables

            // extract the underlying types of the algorithm parameters
            sint tNumSteps = mParameterList.get< sint >( "max_its" );

            // Loop over each variable selected by user - How does the user decide
            // if he/she wants to loop over one variable at a time or perform
            // sweeps over all possible combinations of the selected variables ?
            // In the current setup, the user can perform sweep over only one
            // variable at a time
            for ( uint iv = 0; iv < mNumAdv; ++iv )
            {
                real tOrigValue = mAdvVec(iv,0); // Original variable value

                // Compute the step size of the design variable. Note that this
                // is computed based on the user defined lower and upper bounds.
                real tStepSize = ( mAdvUpVec(iv,0) - mAdvLowVec(iv,0) )/tNumSteps;

                // Loop over number of iterations
                for ( sint k = 0; k <= tNumSteps; ++k )
                {
                    // update adv based upon step size
                    mAdvVec(iv,0) = mAdvLowVec(iv,0) + k*tStepSize;

                    this->func( k );// Perform Forward Analysis
                    this->grad();   // Perform Sensitivity Analysis
                }

                mAdvVec(iv,0) = tOrigValue; // Restore original variable value
            }

            aOptProb = mOptProb; // update aOptProb
        }

        //----------------------------------------------------------------------

        void OptAlgSweep::func( moris::sint aIter )
        {
            // Call to compute objectives and constraints
            OptAlg::func( aIter );
        }


        //----------------------------------------------------------------------
        void OptAlgSweep::grad( )
        {
            // Call to compute derivatives of objectives and constraints
            // w.r.t. advs
            OptAlg::grad( );
        }
    }
}
