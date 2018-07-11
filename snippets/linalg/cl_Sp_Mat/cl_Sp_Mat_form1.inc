moris::Mat< moris::uint > RowInd( 5, 1 );
moris::Mat< moris::uint > ColInd( 5, 1 );
moris::Mat< moris::real > Values( 5, 1 );

RowInd( 0, 0 ) = 2;   RowInd( 1, 0 ) = 1;   RowInd( 2, 0 ) = 2;
RowInd( 3, 0 ) = 3;   RowInd( 4, 0 ) = 0;
ColInd( 0, 0 ) = 0;   ColInd( 1, 0 ) = 1;   ColInd( 2, 0 ) = 1;
ColInd( 3, 0 ) = 2;   ColInd( 4, 0 ) = 3;
Values( 0, 0 ) = 1.0; Values( 1, 0 ) = 2.0; Values( 2, 0 ) = 3.0;
Values( 3, 0 ) = 4.0; Values( 4, 0 ) = 5.0;

moris::Sp_Mat< moris::real > S1( RowInd, ColInd, Values );

