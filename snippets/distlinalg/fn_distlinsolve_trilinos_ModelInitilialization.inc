// ==> assuming 1D bar problem (4 nodes = 3 elements)
moris::uint nd = 4;         		// number of total nodes
moris::uint ne = 3;         		// number of elements
moris::uint ndele = 2;        		// number of nodes per element
moris::uint ndele2= ndele * ndele;  // number of component elemental jacobian


// ==> Trilinos initialization
// Trilinos Parallel Communicator
#ifdef HAVE_MPI
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

// Set the size of the global system to the number of nodes (Global Dimension)
moris::sint NumGlobalElements = nd;

moris::uint index = 0; // Indicate that indexing starts from 0 to the Epetra_Map
// Construct a Map that puts approximately the same number of equations on each processor.
Epetra_Map Map(NumGlobalElements, index, Comm);

// ==>Construct Matricies and vectors for the linear relation Ax=b
// Create a Epetra sparse matrix whose rows have distribution given by the Map.
moris::uint maxRowEntries = 3;
distmat A(Copy, Map, maxRowEntries); // Jacobian Matrix

// Create a Epetra_Vector typedef distvec given by the Map.
distvec b (Map);
distvec x (Map);
