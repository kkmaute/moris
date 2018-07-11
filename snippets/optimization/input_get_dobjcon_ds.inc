void
moris::input::get_dobjcon_ds(
        const moris::Mat< moris::real > &  aAbsDesVarVec,
        const moris::Mat< moris::real > &  aCriteria,
        moris::Mat< moris::real >       & aDObjective_Ds,
        moris::Mat< moris::real >       & aDConstraints_Ds )
{
    // objective = 3*x - x*( crietria = x^2 )
    aDObjective_Ds(0,0) = 3 - aCriteria(0,0);

    // constraint = x - 2
    aDConstraints_Ds(0,0) = 1;
}
