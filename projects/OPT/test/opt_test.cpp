// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "cl_OPT_Interface_Rosenbrock.hpp"
#include "cl_OPT_Problem_Rosenbrock.hpp"
#include "cl_OPT_Manager.hpp"

// ---------------------------------------------------------------------
using namespace moris;

TEST_CASE( "[optimization]" )
{

    SECTION("GCMMA")
    {
        // Create test interface
        moris::opt::Interface_Rosenbrock tInterface;

        // Create problem, assign interface
        moris::opt::Problem_Rosenbrock tProblem(&tInterface);
        tProblem.mConstrained = true;
        tProblem.mFiniteDifferenceObjectives = true;
        tProblem.mFiniteDifferenceConstraints = true;

        // Create algorithm and set parameters
        opt::Algorithm_API tAlgGCMMA("GCMMA");
        tAlgGCMMA.set_param("step_size") = 0.01;    // set the desired step size
        tAlgGCMMA.set_param("penalty")   = 1000.0; // set the desired GCMMA penalty

        // Assign algorithm to cell for manager
        moris::Cell< opt::Algorithm_API > tAlgorithms(1, tAlgGCMMA);

        // Create manager, assign cell of algorithms and problem
        moris::opt::Manager tManager(tAlgorithms, &tProblem);

        // Solve the optimization problem
        tManager.solve_opt_system();

        // Check Solution
        REQUIRE(std::abs(tProblem.get_objectives()(0)) < 1.0e-02); // check value of objective
        REQUIRE(norm(tProblem.get_advs() - 1.0) < 1.0e-02); // check value of design variable, x
    }

    SECTION("SQP")
    {
        // Create test interface
        moris::opt::Interface_Rosenbrock tInterface;

        // Create problem, assign interface
        moris::opt::Problem_Rosenbrock tProblem(&tInterface);
        tProblem.mConstrained = true;

        // Create algorithm and set parameters TODO set parameters for tighter convergence
        opt::Algorithm_API tAlgSQP("SQP");

        // Assign algorithm to cell for manager
        moris::Cell< opt::Algorithm_API > tAlgorithms(1, tAlgSQP);

        // Create manager, assign cell of algorithms and problem
        moris::opt::Manager tManager(tAlgorithms, &tProblem);

        // Solve optimization problem
        tManager.solve_opt_system();

        // Check Solution
        REQUIRE(std::abs(tProblem.get_objectives()(0)) < 1.0e-03); // check value of objective
        REQUIRE(norm(tProblem.get_advs() - 1.0) < 2.0e-02); // check value of design variable, x
    }

//    SECTION("LBFGS")
//    {
//
//        // Create test interface
//        moris::opt::Interface_Rosenbrock tInterface;
//
//        // Create problem, assign interface
//        moris::opt::Problem_Rosenbrock tProblem(&tInterface);
//        tProblem.mConstrained = false;
//
//        // Create algorithm and set parameters
//        opt::Algorithm_API tAlgBFGS("LBFGS");
//
//        // Assign algorithm to cell for manager
//        moris::Cell< opt::Algorithm_API > tAlgorithms(1, tAlgBFGS);
//
//        // Create manager, assign cell of algorithms and problem
//        moris::opt::Manager tManager(tAlgorithms, &tProblem);
//
//        // Solve optimization problem
//        tManager.solve_opt_system();
//
//        moris::print(tProblem.get_advs(), "advs");
//        moris::print(tProblem.get_objectives(), "objective");
//
//
//    }

    SECTION( "SWEEP" )
    {
        // Create test interface
        moris::opt::Interface_Rosenbrock tInterface;

        // Create problem, assign interface
        moris::opt::Problem_Rosenbrock tProblem(&tInterface);
        tProblem.mConstrained = true;

        // Create algorithm and set parameters
        opt::Algorithm_API tAlgSweep("SWEEP");
        tAlgSweep.set_param("num_evaluations_per_adv") = std::string("5, 5");
        std::string tMorisRoot = std::getenv("MORISOUTPUT");
        std::string tHdf5FilePath = tMorisRoot + "sweep.hdf5";
        tAlgSweep.set_param("hdf5_path") = tHdf5FilePath;

        // Assign algorithm to cell for manager
        moris::Cell< opt::Algorithm_API > tAlgorithms(1, tAlgSweep);

        // Create manager, assign cell of algorithms and problem
        moris::opt::Manager tManager(tAlgorithms, &tProblem);

        // Solve optimization problem
        tManager.solve_opt_system();
    }
}

