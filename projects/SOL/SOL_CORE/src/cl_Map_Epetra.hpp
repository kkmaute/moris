/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Map_Epetra.hpp
 *
 */

#pragma once

#include <cstddef>
#include <cassert>
#include <memory>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communicator_Epetra.hpp"
#include "cl_BoostBitset.hpp" // CON/src

#include "cl_SOL_Dist_Map.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_Directory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

class Epetra_MultiVector;

namespace moris
{
    class Solver_Interface;

    class Map_Epetra : public sol::Dist_Map
    {
        private:
            Communicator_Epetra      mEpetraComm;

            Epetra_Map * mEpetraMap             = nullptr;
            Epetra_Map * mEpetraPointMap        = nullptr;
            Epetra_Map * mFullOverlappingMap     = nullptr;   //FIXME move this map to the solver warehouse

            Epetra_MultiVector *mFullToFreePoint = nullptr;

            void build_point_map();

            void translator(
                    const moris::uint      & aNumMaxDofs,
                    const moris::uint      & aNumMyDofs,
                    const Matrix< DDSMat > & aMyLocaltoGlobalMap,
                    Matrix< DDSMat >       & aMyGlobalConstraintDofs,
                    const Matrix< DDUMat > & aMyConstraintDofs );

        protected:

        public:

            //-------------------------------------------------------------------------------------------------------------

            Map_Epetra(
                    const Matrix< DDSMat > & aMyGlobalIds,
                    const Matrix< DDUMat > & aMyConstraintDofs );

            //-------------------------------------------------------------------------------------------------------------

            Map_Epetra( const Matrix< DDSMat > & aMyGlobalIds );

            Map_Epetra( const Vector< sint >& aMyGlobalIds );

            //-------------------------------------------------------------------------------------------------------------
            /** Destructor */
            ~Map_Epetra() override;

            //-------------------------------------------------------------------------------------------------------------

            Epetra_Map * get_epetra_map() override      { return mEpetraMap; };
            Epetra_Map * get_epetra_map() const override{ return mEpetraMap; };

            //-------------------------------------------------------------------------------------------------------------

            Epetra_Map * get_epetra_point_map() override      { return mEpetraPointMap; };
            Epetra_Map * get_epetra_point_map() const override{ return mEpetraPointMap; };

            //-------------------------------------------------------------------------------------------------------------
            moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const override;

            //-------------------------------------------------------------------------------------------------------------

            void build_dof_translator(
                    const Matrix< IdMat > & aFullMap,
                    const bool aFlag ) override;

            //-------------------------------------------------------------------------------------------------------------

            void translate_ids_to_free_point_ids(
                    const moris::Matrix< IdMat > & aIdsIn,
                    moris::Matrix< IdMat >       & aIdsOut,
                    const bool                   & aIsBuildGraph = true ) override;

            void translate_ids_to_free_point_ids(
                    const Vector< sint >& aIdsIn,
                    Vector< sint >&       aIdsOut,
                    bool                  aIsBuildGraph = true ) override;

            //-------------------------------------------------------------------------------------------------------------

            void print() override;

    };
}

