// ==> set solver options for direct solver
KSPCreate(PETSC_COMM_WORLD,&ksp);  // create KSP environment
KSPSetType(ksp,KSPPREONLY);        // build KSP for a particular solver
                                   // KSPPREONLY should be set when direct solver is used.
KSPGetPC(ksp,&pc);                 // return the pointer to the preconditioner
PCSetType(pc,PCLU);                // build preconditioner context, LU preconditioner
KSPSetOperators(ksp,A,A);          // set the matrix associated with the linear system

// ==> direct solver
KSPSolve(ksp,b,x);                 // solve a linear system with direct solver
