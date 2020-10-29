/*
 * cl_SOL_Dist_Map.hpp
 *
 *  Created on: Jun 11, 2018
 *      Author: schmidt
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

#include <petscao.h>

#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace sol
    {
        class Dist_Map
        {
        public:

            // ---------------------------------------------------------------------------------------------------------
            Dist_Map()
            {};

            // ---------------------------------------------------------------------------------------------------------
            /** Destructor */
            virtual ~Dist_Map()
            {};

            // ---------------------------------------------------------------------------------------------------------

            virtual moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const = 0;

            // ---------------------------------------------------------------------------------------------------------

            /**
             * @brief Get Epetra free map.
             *
             * @return  Map object. Either Epetra_Map or AO
             */
            virtual Epetra_Map* get_epetra_map()
            {
                MORIS_ERROR( false, "get_epetra_map() function has no child implementation" );
                return nullptr;
            };

            virtual Epetra_Map* get_epetra_map() const
            {
                MORIS_ERROR( false, "get_epetra_map() function has no child implementation" );
                return nullptr;
            };

            virtual AO get_petsc_map()
            {
                MORIS_ERROR( false, "get_petsc_map() function has no child implementation" );
                return nullptr;
            }
            virtual AO get_petsc_map() const
            {
                MORIS_ERROR( false, "get_petsc_map() function has no child implementation" );
                return nullptr;
            }

        };
    }
}

#endif /* SRC_DISTLINALG_CL_MAP_CLASS_HPP_ */
