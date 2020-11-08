#ifndef SRC_DISTLINALG_CL_MAP_EPETRA_HPP_
#define SRC_DISTLINALG_CL_MAP_EPETRA_HPP_

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

namespace moris
{
    class Map_Epetra : public sol::Dist_Map
    {
        private:
            Communicator_Epetra      mEpetraComm;

            Epetra_Map * mEpetraMap = nullptr;

            void translator( const moris::uint      & aNumMaxDofs,
                    const moris::uint      & aNumMyDofs,
                    const Matrix< DDSMat > & aMyLocaltoGlobalMap,
                    Matrix< DDSMat > & aMyGlobalConstraintDofs,
                    const Matrix< DDUMat > & aMyConstraintDofs );

        protected:

        public:

            //-------------------------------------------------------------------------------------------------------------

            Map_Epetra( const Matrix< DDSMat > & aMyGlobalIds,
                    const Matrix< DDUMat > & aMyConstraintDofs );

            //-------------------------------------------------------------------------------------------------------------

            Map_Epetra( const Matrix< DDSMat > & aMyGlobalIds );

            //-------------------------------------------------------------------------------------------------------------
            /** Destructor */
            ~Map_Epetra();

            //-------------------------------------------------------------------------------------------------------------

            Epetra_Map * get_epetra_map()      { return mEpetraMap; };
            Epetra_Map * get_epetra_map() const{ return mEpetraMap; };

            //-------------------------------------------------------------------------------------------------------------
            moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const;

            //-------------------------------------------------------------------------------------------------------------

            void print();

    };
}

#endif /* SRC_DISTLINALG_CL_MAP_EPETRA_HPP_ */
