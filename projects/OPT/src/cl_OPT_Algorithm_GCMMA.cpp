//Third pary header files
#include "optalggcmmacall.hpp"
#include "mma.hpp"

// MORIS project header files
#include "cl_OPT_Algorithm_GCMMA.hpp" // OPT/src

// -----------------------------------------------------------------------------
using namespace moris;


OptAlgGCMMA::OptAlgGCMMA( ) :
        Algorithm(),
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

void OptAlgGCMMA::solve(moris::opt::Problem* aOptProb )
{
    mProblem = aOptProb; // set the member variable mProblem to aOptProb

    Algorithm::initialize(); // initialize the base class member variables

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
    auto tAdv = mProblem->get_advs().data();
    auto tUpperBounds = mProblem->get_upper_bounds().data();
    auto tLowerBounds = mProblem->get_lower_bounds().data();

    mPrint = false;

    // create an object of type MMAgc solver
    MMAgc mmaAlg(this,
                 tAdv, tUpperBounds, tLowerBounds,
                 mProblem->get_num_advs(), mProblem->get_num_constraints(), tMaxIt, tItsub,
                 tAcc, tSa, tSb, tSc,
                 tDstep, tPenal, NULL, mPrint );

    mResFlag = mmaAlg.solve(); // call the the gcmma solve

    this->printresult(); // print the result of the optimization algorithm

    aOptProb = mProblem; // update aOptProb

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
    // Update the ADV matrix
    Matrix<DDRMat> tADVs(aAdv, aOptAlgGCMMA->mProblem->get_num_advs(), 1);
    aOptAlgGCMMA->mProblem->set_advs(tADVs);

    // Set update for objectives and constraints
    aOptAlgGCMMA->mProblem->mUpdateObjectives = true;
    aOptAlgGCMMA->mProblem->mUpdateConstraints = true;

    // Convert outputs from type MORIS
    aObjval = aOptAlgGCMMA->mProblem->get_objectives()(0);

    // Update the pointer of constraints
    auto tConval = aOptAlgGCMMA->mProblem->get_constraints().data();
    std::copy(tConval, tConval + aOptAlgGCMMA->mProblem->get_num_constraints(), aConval );
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
    aOptAlgGCMMA->mActive = Matrix< DDSMat >  (*aActive, aOptAlgGCMMA->mProblem->get_num_constraints(), 1 );

    // Update the ADV matrix
    Matrix<DDRMat> tADVs(aAdv, aOptAlgGCMMA->mProblem->get_num_advs(), 1);
    aOptAlgGCMMA->mProblem->set_advs(tADVs);

    // Set an update for the gradients
    aOptAlgGCMMA->mProblem->mUpdateObjectiveGradient = true;
    aOptAlgGCMMA->mProblem->mUpdateConstraintGradient = true;

    // Get the objective gradient
    auto tD_Obj = aOptAlgGCMMA->mProblem->get_objective_gradients().data();
    std::copy( tD_Obj, tD_Obj + aOptAlgGCMMA->mProblem->get_num_advs(), aD_Obj );

    // Get the constraint gradient as a MORIS Matrix
    Matrix<DDRMat> tD_Con = aOptAlgGCMMA->mProblem->get_constraint_gradients();

    // Assign to array
    for (moris::uint i = 0; i < aOptAlgGCMMA->mProblem->get_num_constraints(); ++i )
    {
        // loop over number of constraints
        for ( moris::uint j = 0; j < aOptAlgGCMMA->mProblem->get_num_advs(); ++j )
        {
            // Copy data from the matrix to pointer aD_Con
            aD_Con[i][j] = tD_Con(i, j);
        }
    }
}

