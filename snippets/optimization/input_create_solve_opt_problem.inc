void
moris::input::create_solve_opt_problem(
        std::string                 aAlgorithm,
        moris::real &               aObjective,
        moris::Mat< moris::real > & aAdvVec )
{
    // these are user defined variables and depend on the optimization
    // problem the user wishes to solve
    uint aNumMasterAdvs    = 1; // set desired number of advs. To be read from the mesh
    uint aNumCriteria      = 1; // set desired number of criteria
    uint aNumConstraints   = 1; // set total number of constraints
    uint aNumEqConstraints = 0; // set desired number of equality constraints

    // Before we proceed, let's initialize the global variables. NOTE
    // that this is performed specifically to facilitate setting up of
    // unit tests in the current framework and does not constitute a
    // part of how the optimization module operates
    Problem = aProblem;
    ADVs.set_size( aNumMasterAdvs, 1, 0.0 );

    // create an object of type Problem based on the variables defined
    // above
    moris::opt::Problem tOptProb ( aNumMasterAdvs, aNumCriteria, aNumConstraints, aNumEqConstraints );

    tOptProb.mIsOptAnlyt = true; // set the analytical problem flag to true

    // Create an object of type Manager
    opt::Manager aOptManager( tOptProb );

    // solve the optimization problem
    aOptManager.solve_opt_system( );

    // Update the arguments, to be passed to the test file, based on the
    // values of the optimization run
    aObjective = Objective;
    aAdvVec    = ADVs;
}
