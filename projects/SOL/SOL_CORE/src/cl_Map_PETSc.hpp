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
#include "cl_BoostBitset.hpp"            // CON/src
#include "cl_Communication_Tools.hpp"    // COM/src

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
        // flag for full vector maps
        bool mIsFullMap = false;

        // map between moris IDs and Petsc IDs
        AO mPETScMap = nullptr;

        // map from moris IDs to local indices (used by full maps only)
        ISLocalToGlobalMapping mMorisIDtoIndexMap = nullptr;

        // list of Petsc IDs corresponding to moris ids (used by full maps only)
        IS mPETScIDs = nullptr;

        void translator(
                const moris::uint&             aNumMaxDofs,
                const moris::uint&             aNumMyDofs,
                const moris::Matrix< DDSMat >& aMyLocalToGlobalMap,
                moris::Matrix< DDSMat >&       aMyGlobalConstraintDofs,
                const moris::Matrix< DDUMat >& aMyConstraintDofs );

      protected:

      public:
        // constructor for map of owned dofs
        Map_PETSc( const Matrix< DDSMat >& aMyGlobalOwnedIds );

        Map_PETSc( const Vector< sint >& aMyGlobalOwnedIds );

        // constructor for map of owned dofs with constrained dofs
        Map_PETSc(
                const Matrix< DDSMat >& aMyGlobalOwnedIds,
                const Matrix< DDUMat >& aMyConstraintDofs );

        // constructor for map of owned and shared dofs
        Map_PETSc(
                const Matrix< DDSMat >& aMyGlobalOwnedIds,
                const Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds,
                bool                    aIsFullMap );

        // ---------------------------------------------------------------------------------------------------------------

        ~Map_PETSc();

        Matrix< DDSMat >
        map_from_moris_ids_to_indices( const Matrix< DDSMat >& aGlobalIds );

        Vector< sint >
        map_from_moris_ids_to_indices( const Vector< sint >& aGlobalIds );

        // ---------------------------------------------------------------------------------------------------------------

        // check whether map is constructed for full vector
        bool
        is_full_map()
        {
            return mIsFullMap;
        }

        // ---------------------------------------------------------------------------------------------------------------

        // get petsc index list of petsc IDs
        IS
        get_petsc_ids()
        {
            return mPETScIDs;
        }

        // ---------------------------------------------------------------------------------------------------------------

        moris::sint
        return_local_ind_of_global_Id( moris::uint aGlobalId ) const
        {
            MORIS_ERROR( false, "not implemented yet" );

            return -1;
        };

        // ---------------------------------------------------------------------------------------------------------------

        void
        build_dof_translator(
                const Matrix< IdMat >& aFullMap,
                const bool             aFlag )
        {
            MORIS_ERROR( false, "not implemented for petsc yet" );
        };

        // ---------------------------------------------------------------------------------------------------------------

        virtual void
        translate_ids_to_free_point_ids(
                const moris::Matrix< IdMat >& aIdsIn,
                moris::Matrix< IdMat >&       aIdsOut,
                const bool&                   aIsBuildGraph = true )
        {
            MORIS_ERROR( false, "not implemented for petsc yet" );
        };

        void translate_ids_to_free_point_ids(
                const Vector< sint >& aIdsIn,
                Vector< sint >&       aIdsOut,
                bool                  aIsBuildGraph = true )
        {
            MORIS_ERROR( false, "not implemented for petsc yet" );
        }

        // ---------------------------------------------------------------------------------------------------------------

        void
        print()
        {
            AOView( mPETScMap, PETSC_VIEWER_STDOUT_WORLD );
        };

        AO
        get_petsc_map()
        {
            return mPETScMap;
        }

        AO
        get_petsc_map() const
        {
            return mPETScMap;
        }
    };

}    // namespace moris

#endif /* SRC_DISTLINALG_CL_MAP_PETSC_HPP_ */
