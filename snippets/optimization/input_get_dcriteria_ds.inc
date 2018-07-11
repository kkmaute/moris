void
moris::input::get_dcriteria_ds(
        const moris::Mat< moris::real > &  aAbsDesVarVec,
        moris::Mat< moris::real > &        aGradCriteria )
{
    // crietria = x^2
    aGradCriteria(0,0) = 2 * aAbsDesVarVec(0,0);
}
