void
moris::input::initialize(
        moris::Mat< moris::real > &  aAbsDesVarVec,
        moris::Mat< moris::real > &  aAbsDesVarVecUp,
        moris::Mat< moris::real > &  aAbsDesVarVecLow )
{
    aAbsDesVarVec.fill( 1.95 );    // Initialize the advs
    aAbsDesVarVecUp.fill( 2.1 );  // Set upper bounds for the advs
    aAbsDesVarVecLow.fill( 1.9 ); // Set lower bounds for the advs
}
