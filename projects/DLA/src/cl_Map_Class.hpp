/*
 * cl_Map_Class.hpp
 *
 *  Created on: Jun 11, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_MAP_CLASS_HPP_
#define SRC_DISTLINALG_CL_MAP_CLASS_HPP_
// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
 #include <mpi.h>
#endif

#include "linalg.hpp"

#include "Epetra_ConfigDefs.h"
#include "Epetra_Directory.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"

#include <petscao.h>

#include "cl_DistLinAlg_Enums.hpp"

namespace moris
{
class Map_Class
{
private:

protected:

public:
    Epetra_Map * mFreeEpetraMap;
    Epetra_Map * mFullEpetraMap;

    AO          mPETScMap;

// ----------------------------------------------------------------------------------------------------------------------
    Map_Class() :  mFreeEpetraMap(NULL),
                   mFullEpetraMap(NULL)//,
                   //mPETScMap(NULL)
    {
    }

// ----------------------------------------------------------------------------------------------------------------------
     /** Destructor */
    virtual ~Map_Class()
    {
        delete( mFreeEpetraMap );
        delete( mFullEpetraMap );
        //AODestroy( &mPETScMap );
    }
// ----------------------------------------------------------------------------------------------------------------------

    virtual const moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const = 0;

// ----------------------------------------------------------------------------------------------------------------------

    /**
    * @brief Get Epetra Free Map.
    *
    * @return  Map object. Either Epetra_Map or AO
    */
    Epetra_Map* get_epetra_free_map()       {return mFreeEpetraMap;}
    Epetra_Map* get_epetra_free_map() const {return mFreeEpetraMap;}

    /**
     * @brief Get Full Map.
     *
     * @return  Map object. Either Epetra_Map or AO
     */
    Epetra_Map* get_epetra_full_map()       {return mFullEpetraMap;}
    Epetra_Map* get_epetra_full_map() const {return mFullEpetraMap;}

    AO get_petsc_map()       { return mPETScMap; }
    AO get_petsc_map() const { return mPETScMap; }

};
}

#endif /* SRC_DISTLINALG_CL_MAP_CLASS_HPP_ */
