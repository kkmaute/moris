//Third pary header files
#include "optalggcmmacall.hpp"
#include "mma.hpp"

// MORIS project header files
#include "cl_Opt_Alg_GCMMA.hpp" // OPT/src
#include "fn_mem_pointer.hpp" // LNA/src

// -----------------------------------------------------------------------------

OptAlgGCMMA::OptAlgGCMMA( ) :
    OptAlg(),
    mResFlag(0)
{
    // assign default parameter values
    mParameterList.insert( "max_its"      , 100  ); // Allowable GCMMA optimization iterations
    mParameterList.insert( "max_inner_its", 0    ); // Allowable GCMMA inner iterations per every optimization iteration
    mParameterList.insert( "norm_drop"    , 1e-4 ); // GCMMA convergence criteria
    mParameterList.insert( "asymp_adapt0" , 0.5  ); // Initial asymptote adaptation factor
    mParameterList.insert( "asymp_adapt"  , 0.7  ); // Lower asymptote adaptation factor
    mParameterList.insert( "asymp_adaptc" , 1.2  ); // Upper asymptote adaptation factor
    mParameterList.insert( "step_size"    , 0.01 ); // GCMMA step size
    mParameterList.insert( "penalty"      , 100.0); // GCMMA constraint penalty
}

// -----------------------------------------------------------------------------

OptAlgGCMMA::~OptAlgGCMMA()
{
}

// -----------------------------------------------------------------------------

void OptAlgGCMMA::solve( moris::opt::OptProb & aOptProb )
{
    mOptProb = aOptProb; // set the member variable mOptProb to aOptProb

    OptAlg::initialize(); // initialize the base class member variables

    // extract the underlying types of the algorithm parameters and assign
    // to variables that are used to create an object of type MMAgc solver
    moris::sint tMaxIt = mParameterList.get< moris::sint >( "max_its" );
    moris::sint tItsub = mParameterList.get< moris::sint >( "max_inner_its" );
    moris::real tAcc   = mParameterList.get< moris::real >( "norm_drop" );
    moris::real tSa    = mParameterList.get< moris::real >( "asymp_adapt0" );
    moris::real tSb    = mParameterList.get< moris::real >( "asymp_adapt" );
    moris::real tSc    = mParameterList.get< moris::real >( "asymp_adaptc" );
    moris::real tDstep = mParameterList.get< moris::real >( "step_size" );
    moris::real tPenal = mParameterList.get< moris::real >( "penalty" );

    // Note that these pointers are deleted by the the Arma and Eigen
    // libraries themselves.
    auto tAdv    = moris::mem_pointer( mAdvVec );
    auto tAdvUp  = moris::mem_pointer( mAdvUpVec );
    auto tAdvLow = moris::mem_pointer( mAdvLowVec );

    mPrint = false; // default value

    // create an object of type MMAgc solver
    MMAgc mmaAlg( this,
                  tAdv, tAdvUp, tAdvLow,
                  mNumAdv, mNumCon, tMaxIt, tItsub,
                  tAcc, tSa, tSb, tSc,
                  tDstep, tPenal, NULL, mPrint );

    mResFlag = mmaAlg.solve(); // call the the gcmma solve

    printresult(); // print the result of the optimization algorithm

    aOptProb = mOptProb; // update aOptProb

    mmaAlg.cleanup(); // free the memory created by GCMMA
}

// -----------------------------------------------------------------------------

void OptAlgGCMMA::printresult()
{
    if( mPrint )
    {
        std::fprintf( stdout," \nResult of GCMMA\n" );

        switch ( mResFlag )
        {
        case 0:
            std::fprintf( stdout, "\n" );
            std::fprintf( stdout, "THE ALGORITHM HAS CONVERGED.\n" );
            break;
        case 1:
            std::fprintf( stdout, "\n" );
            std::fprintf( stdout, "THE ALGORITHM HAS BEEN STOPPED AFTER MAXIT ITERATIONS.\n" );
            break;
        default:
            std::fprintf( stdout, "\n" );
            std::fprintf( stdout, "Error Message not specified.\n" );
        }
    }
}

// -----------------------------------------------------------------------------

void opt_alg_gcmma_func_wrap(
        OptAlgGCMMA* aOptAlgGCMMA,
        int          aIter,
        double*      aAdv,
        double&      aObjval,
        double*      aConval )
{
    // need to convert to type MORIS to keep coding style consistent
    moris::uint Iter = aIter;

    // Update the ADV matrix
    aOptAlgGCMMA->mAdvVec = moris::Mat< moris::real >( aAdv, aOptAlgGCMMA->mNumAdv, 1 );

    // Call to compute objectives and constraints
    aOptAlgGCMMA->OptAlg::func( Iter );

    // Convert outputs from type MORIS
    aObjval = aOptAlgGCMMA->mObjVal;

    // Update the pointer of constraints
    // Copy data from the matrix mConVal to pointer aConval
    auto tConval = moris::mem_pointer( aOptAlgGCMMA->mConVal );
    std::copy( tConval, tConval + aOptAlgGCMMA->mNumCon, aConval );
}

// -----------------------------------------------------------------------------

void opt_alg_gcmma_grad_wrap(
        OptAlgGCMMA* aOptAlgGCMMA,
        double*      aAdv,
        double*      aD_Obj,
        double**     aD_Con,
        int*         aActive )
{
    // Update the vector of active constraints flag
    aOptAlgGCMMA->mActive = moris::Mat< moris::sint > ( aActive, aOptAlgGCMMA->mNumCon, 1 );

    // Call to compute derivatives of objectives and constraints
    // w.r.t. advs
    aOptAlgGCMMA->OptAlg::grad();

    auto tD_Obj = moris::mem_pointer( aOptAlgGCMMA->mDObj );
    std::copy( tD_Obj, tD_Obj + aOptAlgGCMMA->mNumAdv, aD_Obj );

    // Update the pointer of gradients of objective and constraints
    // loop over number of advs
    for ( moris::uint i=0; i < aOptAlgGCMMA->mNumCon; ++i )
    {
        // loop over number of constraints
        for ( moris::uint j=0; j < aOptAlgGCMMA->mNumAdv; ++j )
        {
            // Copy data from the matrix mDCon to pointer aD_Con
            aD_Con[i][j] = aOptAlgGCMMA->mDCon(i,j);
        }
    }
}

