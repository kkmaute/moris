moris::Mat< moris::cplx > A( 1, 2 );
moris::Mat< moris::cplx > B( 1, 2 );

A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 1.0, 1.0 };

B( 0, 0 ) = { 1.0, 7.0};  B( 0, 1 ) = { 1.0, 1.0 };

moris::Mat< moris::cplx > C( 1, 2 );
C = A % B;