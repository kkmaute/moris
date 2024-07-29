moris::real PI = boost::math::constants::pi<moris::real>();

moris::Mat< moris::real > translation1( { { 0.0, 1.0 } } );
moris::Mat< moris::real > translation2( { { 2.0, 0.0 } } );

moris::Local_Cartesian   CrtSys(translation1);
moris::Local_Cylindrical CylSys(translation2);

moris::Mat< moris::real > testCoords = {
        {0.0, 0.0   },    // glb id 1  // local cart
        {0.0, 1.0   },    // glb id 2  // local cart
        {1.0, 1.0   },    // glb id 7  // local cart
        {1.0, 0.0   },    // glb id 8  // local cart
        {2.0, 0.0   },    // glb id 9  // local cart
        {2.0, 1.0   },    // glb id 10 // local cart
        {2.0, PI/4. },    // glb id 16 // local cyl
        {1.0, PI/4. },    // glb id 17 // local cyl
        {1.0, 0.0   },    // glb id 18 // local cyl
        {2.0, 0.0   } };  // glb id 19 // local cyl

Vector<moris::Loc_Coord_Sys*> testVector =
    { &CrtSys, &CrtSys, &CrtSys, &CrtSys, &CrtSys,
      &CrtSys, &CylSys, &CylSys, &CylSys, &CylSys };

moris::mesh::Mesh_Geometry testGeom(testCoords, testVector );
