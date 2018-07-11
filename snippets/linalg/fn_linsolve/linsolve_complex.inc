moris::Mat< moris::cplx > A( 2, 2 );
moris::Mat< moris::cplx > B( 2, 1 );

A( 0, 0 ) = { 1.0, 1.0 }; A( 0, 1 ) = { 2.0, -1.0 };
A( 1, 0 ) = { 7.0, 0.0 }; A( 1, 1 ) = { 8.0, -2.0 };

B( 0, 0 ) = { 2.0,  7.0 };
B( 1, 0 ) = { 4.0, -9.0 };

moris::Mat< moris::cplx > C = moris::solve( A,B );