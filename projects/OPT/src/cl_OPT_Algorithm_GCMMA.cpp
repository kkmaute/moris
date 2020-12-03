// MORIS
#include "cl_OPT_Algorithm_GCMMA.hpp"
#include "cl_Communication_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

// Third party header files
#include "optalggcmmacall.hpp"
#include "mma.hpp"

using namespace moris;

//----------------------------------------------------------------------------------------------------------------------

OptAlgGCMMA::OptAlgGCMMA(ParameterList aParameterList)
: mRestartIndex(aParameterList.get< moris::sint >( "restart_index" )),
  mMaxIterations(aParameterList.get< moris::sint >( "max_its" )),
  mMaxInnerIterations(aParameterList.get< moris::sint >( "max_inner_its" )),
  mNormDrop(aParameterList.get< moris::real >( "norm_drop" )),
  mAsympAdapt0(aParameterList.get< moris::real >( "asymp_adapt0" )),
  mAsympShrink(aParameterList.get< moris::real >( "asymp_adaptb" )),
  mAsympExpand(aParameterList.get< moris::real >( "asymp_adaptc" )),
  mStepSize(aParameterList.get< moris::real >( "step_size" )),
  mPenalty(aParameterList.get< moris::real >( "penalty" ))
{
}

//----------------------------------------------------------------------------------------------------------------------

OptAlgGCMMA::~OptAlgGCMMA()
{
}

//----------------------------------------------------------------------------------------------------------------------

void OptAlgGCMMA::solve( uint aCurrentOptAlgInd, std::shared_ptr<moris::opt::Problem> aOptProb )
{
    // Trace optimization
    Tracer tTracer( "OptimizationAlgorithm", "GCMMA", "Solve" );

    mCurrentOptAlgInd = aCurrentOptAlgInd;  // set index of current optimization algorithm
    mProblem          = aOptProb;           // set the member variable mProblem to aOptProb

    // Set optimization iteration index for restart
    if ( mRestartIndex > 0)
    {
        gLogger.set_opt_iteration( mRestartIndex );
    }

    // Solve optimization problem
    if (par_rank() == 0)
    {
        // Run gcmma algorithm
        this->gcmma_solve();

        // Communicate that optimization has finished
        mRunning = false;
        this->communicate_running_status();

    }
    else
    {
        // Don't print from these processors
        mPrint = false;

        // Run dummy solve
        this->dummy_solve();
    }

    this->printresult(); // print the result of the optimization algorithm

    aOptProb = mProblem; // update aOptProb
}

//----------------------------------------------------------------------------------------------------------------------

void OptAlgGCMMA::gcmma_solve()
{

    mPrint = false; // FIXME parameter list

    // Note that these pointers are deleted by the the Arma and Eigen
    // libraries themselves.
    auto tAdv         = mProblem->get_advs().data();
    auto tUpperBounds = mProblem->get_upper_bounds().data();
    auto tLowerBounds = mProblem->get_lower_bounds().data();

    // create an object of type MMAgc solver
    MMAgc mmaAlg(this, tAdv, tUpperBounds, tLowerBounds, mProblem->get_num_advs(), mProblem->get_num_constraints(),
            mMaxIterations, mMaxInnerIterations, mNormDrop, mAsympAdapt0, mAsympShrink, mAsympExpand,
            mStepSize, mPenalty, NULL, mPrint);

    mResFlag = mmaAlg.solve(); // call the the gcmma solve

    mmaAlg.cleanup(); // free the memory created by GCMMA
}

//----------------------------------------------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------------------------------------

void opt_alg_gcmma_func_wrap(
        OptAlgGCMMA* aOptAlgGCMMA,
        int          aIter,
        double*      aAdv,
        double&      aObjval,
        double*      aConval )
{
    // Proc 0 needs to communicate that it is still running
    aOptAlgGCMMA->communicate_running_status();

    // Update the ADV matrix
    const Matrix<DDRMat> tADVs(aAdv, aOptAlgGCMMA->mProblem->get_num_advs(), 1);

    // Write restart file
    aOptAlgGCMMA->write_advs_to_file(aIter,tADVs);

    // Recruit help from other procs and solve for criteria
    aOptAlgGCMMA->criteria_solve(tADVs);

    // Set update for objectives and constraints
    aOptAlgGCMMA->mProblem->mUpdateObjectives  = true;
    aOptAlgGCMMA->mProblem->mUpdateConstraints = true;

    // Convert outputs from type MORIS
    aObjval = aOptAlgGCMMA->mProblem->get_objectives()(0);

    // Update the pointer of constraints
    auto tConval = aOptAlgGCMMA->mProblem->get_constraints().data();

    std::copy(tConval, tConval + aOptAlgGCMMA->mProblem->get_num_constraints(), aConval );
}

//----------------------------------------------------------------------------------------------------------------------

void opt_alg_gcmma_grad_wrap(
        OptAlgGCMMA* aOptAlgGCMMA,
        double*      aAdv,
        double*      aD_Obj,
        double**     aD_Con,
        int*         aActive )
{
    // Update the vector of active constraints flag
    aOptAlgGCMMA->mActive = Matrix< DDSMat >(*aActive, aOptAlgGCMMA->mProblem->get_num_constraints(), 1 );

    // Set an update for the gradients
    aOptAlgGCMMA->mProblem->mUpdateObjectiveGradients  = true;
    aOptAlgGCMMA->mProblem->mUpdateConstraintGradients = true;

    // copy objective gradient
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

//----------------------------------------------------------------------------------------------------------------------
