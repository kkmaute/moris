moris::Mat< moris::real > aA2( 3, 3 );
moris::Mat< moris::real > aU;
moris::Mat< moris::real > aS2;
moris::Mat< moris::real > aV;

aA2( 0, 0 ) = 1.0; aA2( 0, 1 ) = 2.0; aA2( 0, 2 ) = 3.0;
aA2( 1, 0 ) = 4.0; aA2( 1, 1 ) = 5.0; aA2( 1, 2 ) = 6.0;
aA2( 2, 0 ) = 7.0; aA2( 2, 1 ) = 8.0; aA2( 2, 2 ) = 9.0;

moris::svd( aU, aS2, aV, aA2 );

moris::Mat< moris::cplx > aAc2( 3, 3 );
moris::Mat< moris::cplx > aUc;
moris::Mat< moris::real > aSc2;
moris::Mat< moris::cplx > aVc;

aAc2( 0, 0 ) = {1.0, 2.0}; aAc2( 0, 1 ) = {3.0, 7.0}; aAc2( 0, 2 ) = {5.0, 9.0};
aAc2( 1, 0 ) = {7.0, 1.0}; aAc2( 1, 1 ) = {9.0, 5.0}; aAc2( 1, 2 ) = {2.0, 8.0};
aAc2( 2, 0 ) = {4.0, 6.0}; aAc2( 2, 1 ) = {6.0, 3.0}; aAc2( 2, 2 ) = {8.0, 1.0};

moris::svd( aUc, aSc2, aVc, aAc2 );
