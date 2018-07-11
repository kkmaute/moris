Cell< opt::OptAlgAPI > define_opt_sol_strategy( )
{
    // create the optimization algorithms
    opt::OptAlgAPI algGCMMA( "GCMMA" );
    opt::OptAlgAPI algSQP  ( "SQP"   );
    opt::OptAlgAPI algSWEEP( "SWEEP" );

    // set the optimization vector of algorithms based on the problem to
    // be solved

    // GCMMA parameters
    algGCMMA.set_param("max_its")   = 10;     // set maximum possible optimization steps
    algGCMMA.set_param("step_size") = 0.001;  // set the desired step size
    algGCMMA.set_param("penalty")   = 1000.0; // set the desired GCMMA penalty

    // SQP parameters
    algSQP.set_param("Major iterations limit") = 10; // set maximum possible optimization steps

    // Return the appropriate vector of algorithms
    return{ algSQP, algGCMMA };
}
