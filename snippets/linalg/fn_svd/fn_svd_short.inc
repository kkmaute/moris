moris::Mat< moris::real > aA( 3, 3 );
moris::Mat< moris::real > aS;

aA( 0, 0 ) = 1.0; aA( 0, 1 ) = 2.0; aA( 0, 2 ) = 3.0;
aA( 1, 0 ) = 4.0; aA( 1, 1 ) = 5.0; aA( 1, 2 ) = 6.0;
aA( 2, 0 ) = 7.0; aA( 2, 1 ) = 8.0; aA( 2, 2 ) = 9.0;

moris::svd( aS, aA );

moris::Mat< moris::cplx > aAc( 3, 3 );
moris::Mat< moris::real > aSc;

aAc( 0, 0 ) = {1.0, 2.0}; aAc( 0, 1 ) = {3.0, 7.0}; aAc( 0, 2 ) = {5.0, 9.0};
aAc( 1, 0 ) = {7.0, 1.0}; aAc( 1, 1 ) = {9.0, 5.0}; aAc( 1, 2 ) = {2.0, 8.0};
aAc( 2, 0 ) = {4.0, 6.0}; aAc( 2, 1 ) = {6.0, 3.0}; aAc( 2, 2 ) = {8.0, 1.0};

moris::svd( aSc, aAc );
