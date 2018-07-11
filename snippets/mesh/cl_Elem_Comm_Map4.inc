// Empty elemental communication map.
moris::Mat<moris::uint> elems3;
moris::Mat<moris::uint> sides3;


moris::mesh::Elem_Comm_Map ecmap3(elems3, sides3);

// Error is produced: elems3 and sides3 have been constructed, but are empty:
moris::size_t num_elems_ecmap3 = ecmap3.num_elems();

