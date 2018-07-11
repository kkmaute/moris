// Amesos stores the timings for different components of the solve in a ParameterList.
ParameterList TimingsList;
Solver->GetTiming (TimingsList);

//
// You can find out how much time was spent in ...
//

// 1) The symbolic factorization
//    (this parameter doesn't always exist; if it doesn't, the
//    two-argument version of get() will return the value of the
//    second argument, which here is zero)
const double sfact_time = TimingsList.get ("Total symbolic factorization time", 0.0);

// 2) The numeric factorization
//    (always exists if NumericFactorization() is called)
const double nfact_time = getParameter<double> (TimingsList, "Total numeric factorization time");

// 3) Solving the linear system
//    (always exists if Solve() is called)
const double solve_time = getParameter<double> (TimingsList, "Total solve time");

// 4) Converting the matrix to the specific solver's input format
//    (always exists if SymbolicFactorization() is called)
const double mtx_conv_time = getParameter<double> (TimingsList, "Total matrix conversion time");

// 5) Redistributing the matrix for each solve to the accepted
//    format for the solver.  (If the matrix to factor is
//    distributed over processors, and the specific solver does not
//    know how to operate on distributed data, then Amesos has to
//    gather the matrix onto one node for the factorization.)
//
//    This may not exist in TimingsList if the matrix didn't need to
//    be redistributed, which is why we use the two-argument form of
//    get() that returns the second argument if the parameter with
//    the given name doesn't exist.
const double mtx_redist_time = TimingsList.get ("Total matrix redistribution time", 0.0);

// 6) Redistributing the vector(s) for each solve to the accepted
//    format for the solver.
//
//    This may not exist in TimingsList if the vector(s) didn't need
//    to be redistributed, which is why we use the two-argument form
//    of get() that returns the second argument if the parameter
//    with the given name doesn't exist.
const double vec_redist_time = TimingsList.get ("Total vector redistribution time", 0.0);

// Just for fun, print out the timings here.  Recall, though, that
// we've already asked Amesos to print out both timings (via the
// "PrintTiming" solver parameter) and "status" (via the
// "PrintStatus" solverparameter).  "Status" in this case means
// interesting details about the factorization.
if (Comm.MyPID() == 0)
{
    cout << endl
    << "Solver timings (in seconds):" << endl
    << "- Symbolic factorization:       " << sfact_time << endl
    << "- Numeric factorization:        " << nfact_time << endl
    << "- Solving the linear system:    " << solve_time << endl
    << "- Converting the matrix:        " << mtx_conv_time << endl
    << "- Redistributing the matrix:    " << mtx_redist_time << endl
    << "- Redistributing the vector(s): " << vec_redist_time << endl
    << endl;
}
