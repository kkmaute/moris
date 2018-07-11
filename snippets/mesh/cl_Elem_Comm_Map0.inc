// Create Element LIDs for elemental communication map.
moris::Mat<moris::size_t> elems1 = {{0,5,1,4,2,3}};

// Create Side LIDs for elemental communication map.
moris::Mat<moris::size_t> sides1 = {{0,1,0,2,3,0}};

// Create element communication map.
moris::mesh::Elem_Comm_Map ecmap1(elems1, sides1);
