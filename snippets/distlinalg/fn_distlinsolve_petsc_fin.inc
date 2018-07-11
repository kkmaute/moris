// ==> deallocate vectors and matrix
VecDestroy(&x);   // destroy Petsc solution vector
VecDestroy(&b);   // destroy Petsc residual vector
MatDestroy(&A);   // destroy Petsc jacobian matrix

// ==> destroy KSP structure
KSPDestroy(&ksp); // destroy KSP environment

// ==> finalization
PetscFinalize();  // finish Petsc solver
