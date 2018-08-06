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

#include "linalg.hpp"
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_Map_Class.hpp"
#include "cl_Solver_Input.hpp"

#include <petscao.h>
#include <petscviewer.h>
#include <petsclog.h>

class Map_PETSc : public moris::Map_Class
{
private:
    void translator( const moris::uint        & aNumMyDofs,
                     const moris::uint        & aNumGlobalDofs,
                     const moris::Mat< int >  & aMyLocaltoGlobalMap,
                           moris::Mat< int >  & aMyGlobalConstraintDofs,
                     const moris::Mat< uint > & aMyConstraintDofs );

protected:

public:
    Map_PETSc( const moris::uint        & aNumMyDofs,
               const moris::Mat<int>    & aMyLocaltoGlobalMap,
               const moris::Mat< uint > & aMyConstraintDofs );

// ---------------------------------------------------------------------------------------------------------------

    ~Map_PETSc();
};


#endif /* SRC_DISTLINALG_CL_MAP_PETSC_HPP_ */
