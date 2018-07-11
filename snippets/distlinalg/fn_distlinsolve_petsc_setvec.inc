// ==> allocation of Petsc vectors
VecCreate(PETSC_COMM_WORLD,&x); // create solution vector
VecSetSizes(x,PETSC_DECIDE,nd); // set vector size (in this case, 4x1 vector)
VecSetFromOptions(x);           // configure the vector from the options database
VecDuplicate(x,&b);             // allocate residual vector by duplicating solution vector

// ==> settings for parallel computing
/*
PETSc parallel vectors are partitioned
by contiguous chunks of rows across the processors.
Determine which vector are locally owned.
*/
//VecGetOwnershipRange(x,&rstart,&rend);
//VecGetLocalSize(x,&nlocal);

// ==> initialization of Petsc solution vector
VecSet(x,zero);                 // set zero vector as the local initial solution vector
VecAssemblyBegin(x);            // assemblye local solution vector into global vector (parallel)
VecAssemblyEnd(x);              // finish assembly part (parallel)
