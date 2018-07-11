// C++ header files.
#include <memory>
#include <iostream>

// MORIS project header files.
#include "Opt_Input.hpp"
#include "typedefs.hpp" // COR/src
#include "cl_Opt_Manager.hpp" // OPT/src
#include "ios.hpp"

namespace moris
{
    namespace opt_input
    {
        uint        Problem;
        real        Objective;
        Mat< real > ADVs;

        // ---------------------------------------------------------------------

        void create_solve_opt_problem(
                uint          aProblem,
                real        & aObjective,
                Mat< real > & aAdvVec )
        {
            // these are user defined variables and depend on the optimization
            // problem the user wishes to solve
            uint aNumMasterAdvs    = 0; // set desired number of advs. To be read from the mesh
            uint aNumCriteria      = 0; // set desired number of criteria
            uint aNumConstraints   = 0; // set total number of constraints
            uint aNumEqConstraints = 0; // set desired number of equality constraints

            // these are user defined variables and depend on the optimization
            // problem the user wishes to solve
            switch ( aProblem )
            {
            case 1:
            case 2:
            case 3:
            case 5:
                aNumMasterAdvs    = 1;
                aNumCriteria      = 1;
                aNumConstraints   = 1;
                break;

            case 4:
                aNumMasterAdvs    = 1;
                aNumCriteria      = 1;
                aNumConstraints   = 0;
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }

            // Before we proceed, let's initialize the global variables. NOTE
            // that this is performed specifically to facilitate setting up of
            // unit tests in the current framework and does not constitute a
            // part of how the optimization module operates
            Problem = aProblem;
            ADVs.set_size( aNumMasterAdvs, 1, 0.0 );

            // create an object of type OptProb based on the variables defined
            // above
            moris::opt::OptProb tOptProb ( aNumMasterAdvs, aNumCriteria, aNumConstraints, aNumEqConstraints );

            tOptProb.mIsOptAnlyt = true; // set the analytical problem flag to true

            // Create an object of type OptManager
            opt::OptManager aOptManager( tOptProb );

            // solve the optimization problem
            aOptManager.solve_opt_system( );

            // Update the arguments, to be passed to the test file, based on the
            // values of the optimization run
            aObjective = Objective;
            aAdvVec    = ADVs;
        }

        // ---------------------------------------------------------------------

        Cell< opt::OptAlgAPI > define_opt_sol_strategy( )
        {
            // create the optimization algorithms
            opt::OptAlgAPI algGCMMA( "GCMMA" );
            opt::OptAlgAPI algSQP  ( "SQP"   );
            opt::OptAlgAPI algSWEEP( "SWEEP" );

            // set the optimization vector of algorithms based on the problem to
            // be solved
            switch ( Problem )
            {
            case 1:
                // GCMMA parameters
                algGCMMA.set_param("max_its")   = 200;    // set maximum possible optimization steps
                algGCMMA.set_param("step_size") = 0.1;    // set the desired step size
                algGCMMA.set_param("penalty")   = 1000.0; // set the desired GCMMA penalty
                break;

            case 2:
                // SQP parameters
                algSQP.set_param("Major iterations limit") = 10; // set maximum possible optimization steps
                break;

            case 3:
                // GCMMA parameters
                algGCMMA.set_param("max_its")   = 10;     // set maximum possible optimization steps
                algGCMMA.set_param("step_size") = 0.001;  // set the desired step size
                algGCMMA.set_param("penalty")   = 1000.0; // set the desired GCMMA penalty

                // SQP parameters
                algSQP.set_param("Major iterations limit") = 10; // set maximum possible optimization steps
                break;

            case 5:
                // SWEEP parameters
                algSWEEP.set_param("max_its") = 300; // set maximum possible optimization steps
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }

            // Return the appropriate vector of algorithms
            if( Problem == 1 )
                return{ algGCMMA };
            else if( Problem == 2 )
                return{ algSQP };
            else if( Problem == 3 )
                return{ algSQP, algGCMMA };
            else
                return{ algSWEEP };
        }

        // ---------------------------------------------------------------------

        void define_advs(
                Mat< real > & aAbsDesVarVec,
                Mat< real > & aAbsDesVarVecUp,
                Mat< real > & aAbsDesVarVecLow )
        {
            // Fill adv related vectors based on problem Id
            switch ( Problem )
            {
            case 1:
            case 2:
            case 3:
                aAbsDesVarVec.fill( 1.95 );   // Initialize the advs
                aAbsDesVarVecUp.fill( 2.1 );  // Set upper bounds for the advs
                aAbsDesVarVecLow.fill( 1.9 ); // Set lower bounds for the advs
                break;

            case 4:
                aAbsDesVarVec.fill( 2.5 );    // Initialize the advs
                aAbsDesVarVecUp.fill( 4.0 );  // Set upper bounds for the advs
                aAbsDesVarVecLow.fill( 1.0 ); // Set lower bounds for the advs
                break;

            case 5:
                aAbsDesVarVec.fill( 2.0 );    // Initialize the advs
                aAbsDesVarVecUp.fill( 2.1 );  // Set upper bounds for the advs
                aAbsDesVarVecLow.fill( 1.9 ); // Set lower bounds for the advs
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        void define_opt_criteria( Mat< real > & aOptCriteria )
        {
            // Specific to this test. Will depend on the problem type
            // Will involving setting up the criteria vector to have knowledge
            // of what forward analysis is to be performed to obtain a
            // particular criteria
        }

        // ---------------------------------------------------------------------

        void get_criteria(
                const Mat< real > & aAbsDesVarVec,
                Mat< real >       & aCriteria )
        {
            switch ( Problem )
            {
            case 1:
            case 2:
            case 3:
            case 5:
                // criteria = x^2
                aCriteria(0,0) = aAbsDesVarVec(0,0) * aAbsDesVarVec(0,0);
                break;

            case 4:
                // criteria = sin(x)
                aCriteria(0,0) = std::sin( aAbsDesVarVec(0,0) );
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        void get_dcriteria_ds(
                const Mat< real > & aAbsDesVarVec,
                Mat< real >       & aGradCriteria )
        {
            switch ( Problem )
            {
            case 1:
            case 2:
            case 3:
            case 5:
                // criteria = x^2
                aGradCriteria(0,0) = 2 * aAbsDesVarVec(0,0);
                break;

            case 4:
                // criteria = sin(x)
                aGradCriteria(0,0) = std::cos( aAbsDesVarVec(0,0) );
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        void get_obj_con(
                const uint          aIter,
                const Mat< real > & aAbsDesVarVec,
                const Mat< real > & aCriteria,
                real              & aObjective,
                Mat< real >       & aConstraints,
                Mat< sint >       & aTypeCon )
        {
            switch ( Problem )
            {
            case 1:
            case 2:
            case 3:
            case 5:
                // objective = 3*x - x*( criteria = x^2 )
                aObjective = 3*aAbsDesVarVec(0,0) - aAbsDesVarVec(0,0)*aCriteria(0,0);

                // constraint = x - 2
                aConstraints(0,0) = aAbsDesVarVec(0,0) - 2;

                // set the type of constraint flags
                if ( aIter == 0 )
                    aTypeCon(0,0) = 1;

                break;

            case 4:
                // objective = 2*( criteria = sin(x) ) - (x^2)/10
                aObjective = 2*aCriteria(0,0) - std::pow( aAbsDesVarVec(0,0), 2 )/10;
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }

            // update the global variables such that when they are passed to the
            // unit test, they have the values from the last optimziation step.
            Objective = aObjective;
            ADVs      = aAbsDesVarVec;
        }

        // ---------------------------------------------------------------------

        void get_dobjcon_ds(
                const Mat< real > & aAbsDesVarVec,
                const Mat< real > & aCriteria,
                Mat< real >       & aDObjective_Ds,
                Mat< real >       & aDConstraints_Ds )
        {
            // compute derivatives for explicit dependency
            switch ( Problem )
            {
            case 1:
            case 2:
            case 3:
            case 5:
                // objective = 3*x - x*( criteria = x^2 )
                aDObjective_Ds(0,0) = 3 - aCriteria(0,0);

                // constraint = x - 2
                aDConstraints_Ds(0,0) = 1;
                break;

            case 4:
                // objective = 2*( criteria = sin(x) ) - (x^2)/10
                aDObjective_Ds(0,0) = -2*aAbsDesVarVec(0,0)/10;
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        void get_dobjcon_dcrit(
                Mat< real >   aAbsDesVarVec,
                Mat< real >   aCriteria,
                Mat< real > & aDObjective_DCrit,
                Mat< real > & aDConstraints_DCrit )
        {
            switch ( Problem )
            {
            case 1:
            case 2:
            case 3:
            case 5:
                // objective = 3*x - x*( criteria = x^2 )
                aDObjective_DCrit (0,0) = -aAbsDesVarVec(0,0);
                break;

            case 4:
                aDObjective_DCrit (0,0) = 2.0;
                break;

            default:
                MORIS_LOG_ERROR << "Problem not yet implemented.";
                assert::error( "In opt_input.cpp" );
            }
        }

        // ---------------------------------------------------------------------
    }
}

