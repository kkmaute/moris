/*
 * cl_MSI_Solver_Interface.cpp
 *
 *  Created on: Sep 23, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_SOL_Dist_Vector.hpp"

namespace moris
{
    namespace MSI
    {

//-------------------------------------------------------------------------------------------------------

    void MSI_Solver_Interface::get_exact_solution_from_hdf5_and_calculate_error( const char* aFilename )
    {
        mPrevSolutionVector->read_vector_from_HDF5( aFilename );

        mSolutionVector->vec_plus_vec( 1.0, *mPrevSolutionVector, -1.0 );
    }

//-------------------------------------------------------------------------------------------------------

    void MSI_Solver_Interface::get_residual_vector_for_output( const char* aFilename )
    {
//        mPrevSolutionVector->read_vector_from_HDF5( aFilename );
        mPrevSolutionVector->vec_put_scalar( 5.0 );

        mSolutionVector->vec_plus_vec( 1.0, *mPrevSolutionVector, 0.0 );
    }

//-------------------------------------------------------------------------------------------------------

    void MSI_Solver_Interface::write_solution_to_hdf5_file( const char* aFilename )
    {
        mSolutionVector->save_vector_to_HDF5( aFilename );
    }

//-------------------------------------------------------------------------------------------------------

    void MSI_Solver_Interface::initiate_output( const uint aOutputIndex,
                                                const real aTime,
                                                const bool aEndOfTimeIteration )
    {
        // end of time iteration that the exodus file should be closed
        mModel->output_solution( aOutputIndex, aTime, aEndOfTimeIteration );
    }

//------------------------------------------------------------------------------

    void MSI_Solver_Interface::set_solution_vector( Dist_Vector * aSolutionVector )
    {
        mSolutionVector = aSolutionVector;
        mMSI->mEquationModel->set_solution_vector( mSolutionVector );
    }

//------------------------------------------------------------------------------

    void MSI_Solver_Interface::set_solution_vector_prev_time_step( Dist_Vector * aSolutionVector )
    {
        mPrevSolutionVector = aSolutionVector;
        mMSI->mEquationModel->set_previous_solution_vector( mPrevSolutionVector );
    }

//------------------------------------------------------------------------------

     void MSI_Solver_Interface::set_adjoint_solution_vector( Dist_Vector * aSolutionVector )
     {
         mSensitivitySolutionVector = aSolutionVector;
         mMSI->mEquationModel->set_adjoint_solution_vector( mSensitivitySolutionVector );
     }

//------------------------------------------------------------------------------

    void MSI_Solver_Interface::set_time( const Matrix< DDRMat> & aTime )
    {
        mTime = aTime;
        mMSI->mEquationModel->set_time( mTime );
    }

//------------------------------------------------------------------------------

    void MSI_Solver_Interface::set_previous_time( const Matrix< DDRMat> & aTime )
    {
        mPrevTime = aTime;
        mMSI->mEquationModel->set_previous_time( mPrevTime );
    }

//------------------------------------------------------------------------------

    const moris::Cell < moris::Matrix< DDRMat> > & MSI_Solver_Interface::get_criteria( const moris::uint & aMySetInd )
    {
        return mMSI->get_equation_set( aMySetInd )->get_QI();
    }

//------------------------------------------------------------------------------

    void MSI_Solver_Interface::set_requested_IQI_names( const moris::Cell< std::string > & aIQINames )
    {
        mMSI->get_equation_model()->set_requested_IQI_names( aIQINames );
    }

//-------------------------------------------------------------------------------------------------------

    moris::uint MSI_Solver_Interface::get_num_rhs()
    {
        return mMSI->get_equation_model()->get_num_rhs();
    }

//-------------------------------------------------------------------------------------------------------

    }
}
