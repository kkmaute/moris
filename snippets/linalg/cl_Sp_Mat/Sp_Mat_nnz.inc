moris::Sp_Mat< moris::real > A(5,5);

A(1,2) = 1.0;
A(2,3) = 2.0;
A(3,4) = 3.0;
A(4,0) = 4.0;

moris::uint NumNonZero = A.get_nnz();