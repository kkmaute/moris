// The linear problem for Amesos to solve.
Epetra_LinearProblem Problem (&A, &x, &b);

// Trilinos Solver parameters used for catching errors and timing checks
using Teuchos::getParameter;
using Teuchos::ParameterList;
using std::cerr;
using std::cout;
using std::endl;

// Amesos_BaseSolver is the interface for all Amesos solvers. It is a pure virtual class.
// Hence, objects of this class cannot be allocated, and can exist only as pointers or references.
Amesos_BaseSolver* Solver = NULL;

// Initialize the Factory which will be used to create objects that implement the Amesos_BaseSolver interface.
Amesos Factory;

// Specifies the solver. String ``SolverType'' can assume one of the following values:
// - "Lapack"- "Klu" - "Umfpack" - "Pardiso" - "Taucs" - "Superlu" - "Superludist" - "Mumps" - "Dscpack"
std::string SolverType = "Klu";
Solver = Factory.Create (SolverType, Problem);

// Factory.Create() returns NULL if the requested solver is not available.
if (Solver == NULL) {
cerr << "Specified solver \"" << SolverType << "\" is not available." << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
}

// Parameters in the list are set using List.set("parameter-name",ParameterValue);
ParameterList List;
List.set ("PrintTiming", false);	// Output Solver timing information
List.set ("PrintStatus", false);	// Output Solver status information
Solver->SetParameters (List);

// To solve Call: SymbolicFactorization(), then NumericFactorization(), and finally Solve().
if (Comm.MyPID() == 0)
Solver->SymbolicFactorization();

// You can change the matrix values at any point before calling Solver->NumericFactorization().
if (Comm.MyPID() == 0)
Solver->NumericFactorization();

// You can change the solution vector LHS and the right-hand side vector RHS at any time before calling Solver->Solve().
if (Comm.MyPID() == 0)
Solver->Solve(); // Solve the system (Ax=b) for x
