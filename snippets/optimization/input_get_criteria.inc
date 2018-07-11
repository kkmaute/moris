void
moris::input::get_criteria(
        const moris::Mat< moris::real > &  aAbsDesVarVec,
        moris::Mat< moris::real > &        aCriteria )
{
    // crietria = x^2
    aCriteria(0,0) = aAbsDesVarVec(0,0) * aAbsDesVarVec(0,0);
}
