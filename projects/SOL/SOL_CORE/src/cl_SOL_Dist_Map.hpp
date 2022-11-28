/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Dist_Map.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_MAP_CLASS_HPP_
#define SRC_DISTLINALG_CL_MAP_CLASS_HPP_

// MORIS header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_Directory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

#ifdef MORIS_HAVE_PETSC
#include <petscao.h>
#endif

#include "cl_SOL_Enums.hpp"

namespace moris
{
    class Solver_Interface;

    namespace sol
    {
        class Dist_Map
        {
          public:
            // ---------------------------------------------------------------------------------------------------------
            Dist_Map(){};

            // ---------------------------------------------------------------------------------------------------------
            /** Destructor */
            virtual ~Dist_Map(){};

            // ---------------------------------------------------------------------------------------------------------

            virtual moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const = 0;

            // ---------------------------------------------------------------------------------------------------------

            virtual void build_dof_translator(
                    const Matrix< IdMat >& aFullMap,
                    const bool             aFlag ) = 0;

            // ---------------------------------------------------------------------------------------------------------

            virtual void translate_ids_to_free_point_ids(
                    const moris::Matrix< IdMat >& aIdsIn,
                    moris::Matrix< IdMat >&       aIdsOut,
                    const bool&                   aIsBuildGraph = true ) = 0;

            // ---------------------------------------------------------------------------------------------------------

            virtual void print() = 0;

            // ---------------------------------------------------------------------------------------------------------

            /**
             * @brief Get Epetra free map.
             *
             * @return  Map object. Either Epetra_Map or AO
             */
            virtual Epetra_Map*
            get_epetra_map()
            {
                MORIS_ERROR( false, "get_epetra_map() function has no child implementation" );
                return nullptr;
            };

            virtual Epetra_Map*
            get_epetra_map() const
            {
                MORIS_ERROR( false, "get_epetra_map() function has no child implementation" );
                return nullptr;
            };

            virtual Epetra_Map*
            get_epetra_point_map()
            {
                MORIS_ERROR( false, "get_epetra_point_map() function has no child implementation" );
                return nullptr;
            };
            virtual Epetra_Map*
            get_epetra_point_map() const
            {
                MORIS_ERROR( false, "get_epetra_point_map() function has no child implementation" );
                return nullptr;
            };

#ifdef MORIS_HAVE_PETSC
            virtual AO
            get_petsc_map()
            {
                MORIS_ERROR( false, "get_petsc_map() function has no child implementation" );
                return nullptr;
            }
            virtual AO
            get_petsc_map() const
            {
                MORIS_ERROR( false, "get_petsc_map() function has no child implementation" );
                return nullptr;
            }
#endif
        };
    }    // namespace sol
}    // namespace moris

#endif /* SRC_DISTLINALG_CL_MAP_CLASS_HPP_ */
