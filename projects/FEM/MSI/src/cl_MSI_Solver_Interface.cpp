/*
 * cl_MSI_Solver_Interface.cpp
 *
 *  Created on: Sep 23, 2018
 *      Author: schmidt
 */
#include "cl_MSI_Solver_Interface.hpp"

#include "cl_Vector.hpp"

namespace moris
{
    namespace MSI
    {

    void MSI_Solver_Interface::get_exact_solution_from_hdf5_and_calculate_error( const char* aFilename )
    {
        mPrevSolutionVector->read_vector_from_HDF5( aFilename );

        mSolutionVector->vec_plus_vec( 1.0, *mPrevSolutionVector, -1.0 );
    }

    void MSI_Solver_Interface::write_solution_to_hdf5_file( const char* aFilename )
    {
        mSolutionVector->save_vector_to_HDF5( aFilename );
    }

    }
}
