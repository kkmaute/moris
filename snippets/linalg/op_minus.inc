moris::Mat< moris::real > Am( 3, 3 );
moris::Mat< moris::real > Bm( 3, 3 );
moris::Mat< moris::real > Cm( 3, 3 );
moris::Mat< moris::real > Dm( 3, 3 );

Am( 0, 0 ) = 1.0; Am( 0, 1 ) = 1.0; Am( 0, 2 ) = 1.0;
Am( 1, 0 ) = 1.0; Am( 1, 1 ) = 1.0; Am( 1, 2 ) = 1.0;
Am( 2, 0 ) = 1.0; Am( 2, 1 ) = 1.0; Am( 2, 2 ) = 1.0;

Bm( 0, 0 ) = 2.0; Bm( 0, 1 ) = 2.0; Bm( 0, 2 ) = 2.0;
Bm( 1, 0 ) = 2.0; Bm( 1, 1 ) = 2.0; Bm( 1, 2 ) = 2.0;
Bm( 2, 0 ) = 2.0; Bm( 2, 1 ) = 2.0; Bm( 2, 2 ) = 2.0;

Cm( 0, 0 ) = 3.0; Cm( 0, 1 ) = 3.0; Cm( 0, 2 ) = 3.0;
Cm( 1, 0 ) = 3.0; Cm( 1, 1 ) = 3.0; Cm( 1, 2 ) = 3.0;
Cm( 2, 0 ) = 3.0; Cm( 2, 1 ) = 3.0; Cm( 2, 2 ) = 3.0;

Dm = Am - Bm - Cm;

moris::Sp_Mat< moris::real > As( 3, 3 );
moris::Sp_Mat< moris::real > Bs( 3, 3 );
moris::Sp_Mat< moris::real > Cs( 3, 3 );
moris::Sp_Mat< moris::real > Ds( 3, 3 );

As( 0, 0 ) = 1.0; As( 0, 1 ) = 1.0; As( 0, 2 ) = 1.0;
As( 1, 0 ) = 1.0; As( 1, 1 ) = 1.0; As( 1, 2 ) = 1.0;
As( 2, 0 ) = 1.0; As( 2, 1 ) = 1.0; As( 2, 2 ) = 1.0;

Bs( 0, 0 ) = 2.0; Bs( 0, 1 ) = 2.0; Bs( 0, 2 ) = 2.0;
Bs( 1, 0 ) = 2.0; Bs( 1, 1 ) = 2.0; Bs( 1, 2 ) = 2.0;
Bs( 2, 0 ) = 2.0; Bs( 2, 1 ) = 2.0; Bs( 2, 2 ) = 2.0;

Cs( 0, 0 ) = 3.0; Cs( 0, 1 ) = 3.0; Cs( 0, 2 ) = 3.0;
Cs( 1, 0 ) = 3.0; Cs( 1, 1 ) = 3.0; Cs( 1, 2 ) = 3.0;
Cs( 2, 0 ) = 3.0; Cs( 2, 1 ) = 3.0; Cs( 2, 2 ) = 3.0;

Ds = As - Bs - Cs;

moris::Mat< moris::real > Ac( 3, 1 );
moris::Mat< moris::real > Bc( 3, 1 );
moris::Mat< moris::real > Cc( 3, 1 );
moris::Mat< moris::real > Dc( 3, 1 );

Ac( 0, 0 ) = 1.0;
Ac( 1, 0 ) = 1.0;
Ac( 2, 0 ) = 1.0;

Bc( 0, 0 ) = 2.0;
Bc( 1, 0 ) = 2.0;
Bc( 2, 0 ) = 2.0;

Cc( 0, 0 ) = 3.0;
Cc( 1, 0 ) = 3.0;
Cc( 2, 0 ) = 3.0;

Dc = Ac - Bc - Cc;
