#ifndef SRC_DISTLINALG_CL_MAP_EPETRA_HPP_
#define SRC_DISTLINALG_CL_MAP_EPETRA_HPP_

#ifdef MORIS_HAVE_PARALLEL
 #include "Epetra_MpiComm.h"
 #include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include <cstddef>
#include <cassert>
#include <memory>

#include "linalg.hpp"
#include "cl_Communicator_Epetra.hpp"
#include "cl_BoostBitset.hpp" // CON/src

#include "cl_Map_Class.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_Directory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

namespace moris
{
class Map_Epetra : public Map_Class
{
private:
    Communicator_Epetra      mEpetraComm;

    void translator( const moris::uint        & aNumMyDofs,
                     const moris::uint        & aNumGlobalDofs,
                     const moris::Mat< int >  & aMyLocaltoGlobalMap,
                           moris::Mat< int >  & aMyGlobalConstraintDofs,
                     const moris::Mat< uint > & aMyConstraintDofs );

protected:

public:
    Map_Epetra( const moris::uint        & aNumMyDofs,
                const moris::Mat< int >  & aMyLocaltoGlobalMap,
                const moris::Mat< uint > & aMyConstraintDofs );

//-------------------------------------------------------------------------------------------------------------
/** Destructor */
~Map_Epetra();

};
}

#endif /* SRC_DISTLINALG_CL_MAP_EPETRA_HPP_ */
