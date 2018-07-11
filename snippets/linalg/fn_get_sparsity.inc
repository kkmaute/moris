moris::Sp_Mat< moris::real > aMat(5,5);

aMat(1,2) = 1.0;
aMat(2,3) = 2.0;
aMat(3,4) = 3.0;
aMat(4,0) = 4.0;

moris::Mat< moris::uint > Sparsity( aMat.get_nnz(), 2 );

moris::get_sparsity( Sparsity, aMat);