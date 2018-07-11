// ==> allocation of Petsc matrix
MatCreateSeqAIJ(MPI_COMM_SELF,nd,nd,nd,0,&A); // allocate A matrix
