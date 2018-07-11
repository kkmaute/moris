// ==> initialization of Petsc
distvec	x;                    // solution vector (Vec in standard Petsc type)
distvec b;                    // residual vector (Vec in standard Petsc type)
distmat	A;                    // jacobian matrix (Mat in standard Petsc type)
KSP     ksp;                  // linear solver context of KSP
PC      pc;                   // preconditioner context
distint	n,rstart,rend,nlocal; // Petsc integer for parallel computing
distreal zero = 0;            // zero value as Petsc scalar for initialization of x

static char help[] = "Appends to an ASCII file.\n\n";
PetscInitialize((int*)0,(char)0,(char*)0,help); // initialize Petsc
