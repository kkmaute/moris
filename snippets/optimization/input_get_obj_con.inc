void
moris::input::get_obj_con(
        const moris::Mat< moris::real > &  aAbsDesVarVec,
        const moris::Mat< moris::real > &  aCriteria,
        moris::real &                      aObjective,
        moris::Mat< moris::real > &        aConstraints )
{
    // objective = 3*x - x*( crietria = x^2 )
    aObjective = 3*aAbsDesVarVec(0,0) - aAbsDesVarVec(0,0)*aCriteria(0,0);

    // constraint = x - 2
    aConstraints(0,0) = aAbsDesVarVec(0,0) - 2;

    // update the global variables such that when they are passed to the
    // unit test, they have the values from the last optimziation step.
    moris::input::Objective = aObjective;
    moris::input::ADVs      = aAbsDesVarVec;
}