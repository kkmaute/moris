/*
 * MapPETSc.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_MAP_PETSC_HPP_
#define SRC_DISTLINALG_CL_MAP_PETSC_HPP_

#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include <cstddef>
#include <cassert>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_Map_Class.hpp"
#include "cl_Solver_Input.hpp"

#include <petscao.h>
#include <petscviewer.h>
#include <petsclog.h>
namespace moris
{

class Map_PETSc : public moris::Map_Class
{
private:
    void translator( const moris::uint      & aNumMyDofs,
                     const moris::uint      & aNumGlobalDofs,
                     const moris::Matrix< DDSMat > & aMyLocaltoGlobalMap,
                           moris::Matrix< DDSMat > & aMyGlobalConstraintDofs,
                     const moris::Matrix< DDUMat > & aMyConstraintDofs );

protected:

public:
    Map_PETSc( const moris::uint      & aNumMyDofs,
               const Matrix< DDSMat > & aMyLocaltoGlobalMap,
               const Matrix< DDUMat > & aMyConstraintDofs );

// ---------------------------------------------------------------------------------------------------------------

    ~Map_PETSc();

    // ---------------------------------------------------------------------------------------------------------------
    const moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const
    {
        MORIS_ERROR( false, "not implemented yet");

        return -1;
    };
};

}


#endif /* SRC_DISTLINALG_CL_MAP_PETSC_HPP_ */
