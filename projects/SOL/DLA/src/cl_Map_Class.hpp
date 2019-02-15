/*
 * cl_Map_Class.hpp
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

#include "cl_DLA_Enums.hpp"

namespace moris
{
class Map_Class
{
private:

protected:

public:
    Epetra_Map * mFreeEpetraMap;
    Epetra_Map * mFullEpetraMap;
    Epetra_Map * mFullOverlappingEpetraMap;

    AO          mPETScMap;

// ----------------------------------------------------------------------------------------------------------------------
    Map_Class() :  mFreeEpetraMap(NULL),
                   mFullEpetraMap(NULL),
                   mFullOverlappingEpetraMap(NULL),
                   mPETScMap(NULL)
    {
    }

// ----------------------------------------------------------------------------------------------------------------------
     /** Destructor */
    virtual ~Map_Class()
    {
        delete( mFreeEpetraMap );
        delete( mFullEpetraMap );
        delete( mFullOverlappingEpetraMap );
        AODestroy( &mPETScMap );
    }
// ----------------------------------------------------------------------------------------------------------------------

    virtual moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const = 0;

// ----------------------------------------------------------------------------------------------------------------------

    /**
    * @brief Get Epetra free map.
    *
    * @return  Map object. Either Epetra_Map or AO
    */
    Epetra_Map* get_epetra_free_map()       {return mFreeEpetraMap;}
    Epetra_Map* get_epetra_free_map() const {return mFreeEpetraMap;}

    /**
     * @brief Get full map.
     *
     * @return  Map object. Either Epetra_Map or AO
     */
    Epetra_Map* get_epetra_full_map()       {return mFullEpetraMap;}
    Epetra_Map* get_epetra_full_map() const {return mFullEpetraMap;}

    /**
     * @brief Get full overlapping map.
     *
     * @return  Map object. Either Epetra_Map or AO
     */
    Epetra_Map* get_epetra_full_overlapping_map()       {return mFullOverlappingEpetraMap;}
    Epetra_Map* get_epetra_full_overlapping_map() const {return mFullOverlappingEpetraMap;}


    AO get_petsc_map()       { return mPETScMap; }
    AO get_petsc_map() const { return mPETScMap; }

};
}

#endif /* SRC_DISTLINALG_CL_MAP_CLASS_HPP_ */
