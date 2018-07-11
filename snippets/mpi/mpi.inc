// MPI communicator.
boost::mpi::communicator mpi_comm;

// Rank of processor in communicator.
moris::uint mpi_rank = mpi_comm.rank();

// Number of processors in communicator.
moris::uint mpi_size = mpi_comm.size();

