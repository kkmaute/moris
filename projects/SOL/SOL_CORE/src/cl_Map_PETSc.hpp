/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map_PETSc.hpp
 *
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
#include "cl_Communication_Tools.hpp" // COM/src

#include "cl_SOL_Dist_Map.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include <petscao.h>
#include <petscviewer.h>
#include <petsclog.h>
namespace moris
{

class Map_PETSc : public moris::sol::Dist_Map
{
private:

    AO          mPETScMap = nullptr;

    void translator( const moris::uint             & aNumMaxDofs,
                     const moris::uint             & aNumMyDofs,
                     const moris::Matrix< DDSMat > & aMyLocaltoGlobalMap,
                           moris::Matrix< DDSMat > & aMyGlobalConstraintDofs,
                     const moris::Matrix< DDUMat > & aMyConstraintDofs );

protected:

public:
    Map_PETSc( const Matrix< DDSMat > & aMyGlobalIds,
               const Matrix< DDUMat > & aMyConstraintDofs );

    Map_PETSc( const Matrix< DDSMat > & aMyGlobalIds );

// ---------------------------------------------------------------------------------------------------------------

    ~Map_PETSc();

    // ---------------------------------------------------------------------------------------------------------------
    moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const
    {
        MORIS_ERROR( false, "not implemented yet");

        return -1;
    };

    // ---------------------------------------------------------------------------------------------------------------

    void build_dof_translator(
            const Matrix< IdMat > & aFullMap,
            const bool aFlag )
    {
        MORIS_ERROR( false, "not implemented for petsc yet");
    };

    // ---------------------------------------------------------------------------------------------------------------

    virtual void translate_ids_to_free_point_ids(
            const moris::Matrix< IdMat > & aIdsIn,
            moris::Matrix< IdMat >       & aIdsOut,
            const bool                   & aIsBuildGraph = true )
    {
        MORIS_ERROR( false, "not implemented for petsc yet");
    };

    // ---------------------------------------------------------------------------------------------------------------
    void print()
    {
        MORIS_ERROR( false, "not implemented yet");
    };

    AO get_petsc_map()       { return mPETScMap; }
    AO get_petsc_map() const { return mPETScMap; }
};

}

#endif /* SRC_DISTLINALG_CL_MAP_PETSC_HPP_ */

