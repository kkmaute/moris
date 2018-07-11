moris::Mat< moris::real > Am( 3, 3 );
moris::Mat< moris::real > Bm( 3, 3 );
moris::Mat< moris::real > Cm( 3, 3 );

Am( 0, 0 ) = 1.0; Am( 0, 1 ) = 1.0; Am( 0, 2 ) = 1.0;
Am( 1, 0 ) = 1.0; Am( 1, 1 ) = 1.0; Am( 1, 2 ) = 1.0;
Am( 2, 0 ) = 1.0; Am( 2, 1 ) = 1.0; Am( 2, 2 ) = 1.0;

Bm( 0, 0 ) = 2.0; Bm( 0, 1 ) = 2.0; Bm( 0, 2 ) = 2.0;
Bm( 1, 0 ) = 2.0; Bm( 1, 1 ) = 2.0; Bm( 1, 2 ) = 2.0;
Bm( 2, 0 ) = 2.0; Bm( 2, 1 ) = 2.0; Bm( 2, 2 ) = 2.0;

Cm = Am * Bm;

moris::Sp_Mat< moris::real > As( 3, 3 );
moris::Sp_Mat< moris::real > Bs( 3, 3 );
moris::Sp_Mat< moris::real > Cs( 3, 3 );

As( 0, 0 ) = 1.0; As( 0, 1 ) = 1.0; As( 0, 2 ) = 1.0;
As( 1, 0 ) = 1.0; As( 1, 1 ) = 1.0; As( 1, 2 ) = 1.0;
As( 2, 0 ) = 1.0; As( 2, 1 ) = 1.0; As( 2, 2 ) = 1.0;

Bs( 0, 0 ) = 2.0; Bs( 0, 1 ) = 2.0; Bs( 0, 2 ) = 2.0;
Bs( 1, 0 ) = 2.0; Bs( 1, 1 ) = 2.0; Bs( 1, 2 ) = 2.0;
Bs( 2, 0 ) = 2.0; Bs( 2, 1 ) = 2.0; Bs( 2, 2 ) = 2.0;

Cs = As * Bs;