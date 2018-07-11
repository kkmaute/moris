void
moris::input::get_dobjcon_dcrit(
        moris::Mat< moris::real >   aAbsDesVarVec,
        moris::Mat< moris::real >   aCriteria,
        moris::Mat< moris::real > & aDObjective_DCrit,
        moris::Mat< moris::real > & aDConstraints_DCrit)
{
    // objective = 3*x - x*( criteria = x^2 )
    aDObjective_DCrit (0,0) = -aAbsDesVarVec(0,0);
}
