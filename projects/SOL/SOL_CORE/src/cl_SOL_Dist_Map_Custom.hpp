/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SOL_Dist_Map_Custom.hpp
 *
 */

#pragma once

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

#include "cl_Communicator_Epetra.hpp"
#include "Epetra_Map.h"

namespace moris
{
    class Dist_Map_Custom : public sol::Dist_Map
    {
    public:
        // this is a custom AO object in petsc
        std::unordered_map< sint, sint > mApplicationToPetsc;

        // flag for full vector maps
        bool mIsFullMap = false;

        // map from moris IDs to local indices (used by full maps only)
        ISLocalToGlobalMapping mMorisIDtoIndexMap = nullptr;

        // list of owned and shared moris IDs ( both for full and not full maps)
        Matrix<DDSMat>  mMorisIDsOwnedAndShared;

        // list of Petsc IDs corresponding to moris ids (used by full maps only)
        Matrix<DDSMat> mMorisIDsOwned;
        
        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        //constrcutor
        Dist_Map_Custom(const Matrix< DDSMat >& aMyGlobalOwnedIds,const Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds);

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        //constrcutor
        Dist_Map_Custom(const Matrix< DDSMat >& aMyGlobalOwnedIds,const Matrix< DDSMat >& aMyGlobalOwnedAndSharedIds, bool );

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        Dist_Map_Custom(const Matrix< DDSMat >& aMyGlobalOwnedIds,const Matrix< DDUMat >& aMyGlobalOwnedAndSharedIds); 

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
         Dist_Map_Custom(const Matrix< DDSMat >& aMyGlobalOwnedIds){
            MORIS_ERROR( false, "Dist_Map_Custom::Dist_Map_Custom - This constructor is not implemented" );
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             Dist_Map_Custom(const Vector< sint >& aMyGlobalOwnedIds){
            MORIS_ERROR( false, "Dist_Map_Custom::Dist_Map_Custom - This constructor is not implemented" );
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        //destructor
        ~Dist_Map_Custom() override{}

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // retuern the local index of the global id
     void
    map_from_moris_ids_to_petsc_ids( Matrix< DDSMat >& aGlobalIds );

    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
         void
    map_from_moris_ids_to_petsc_ids( Matrix< DDUMat >& aGlobalIds );

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        bool
        is_full_map()
        {
            return mIsFullMap;
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        void
        print() override
        {
            // since all the processors have the same data, proc 0 prints it
            if ( par_rank() == 0 )
            {
                std::cout << " moris id" << " => " << "petsc id" << '\n';
                for(auto it = mApplicationToPetsc.begin(); it != mApplicationToPetsc.end(); ++it)
                {
                    std::cout << it->first << " => " << it->second << '\n';
                }
            }
            barrier();
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        Matrix< DDSMat > const &
        get_moris_ids_owned_and_shared() const
        {
            return mMorisIDsOwnedAndShared;
        }

        //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        Matrix< DDSMat > const &
        get_moris_ids_owned() const
        {
            return mMorisIDsOwned;
        }

                   // ---------------------------------------------------------------------------------------------------------

            moris::sint return_local_ind_of_global_Id( moris::uint aGlobalId ) const override
            {
                MORIS_ERROR( false, "return_local_ind_of_global_Id() function has no child implementation" );
                return -1;
            }

            // ---------------------------------------------------------------------------------------------------------

            void build_dof_translator(
                    const Matrix< IdMat >& aFullMap,
                    const bool             aFlag ) override
                    {
                        MORIS_ERROR( false, "build_dof_translator() function has no child implementation" );
                    }

            // ---------------------------------------------------------------------------------------------------------

            void translate_ids_to_free_point_ids(
                    const moris::Matrix< IdMat >& aIdsIn,
                    moris::Matrix< IdMat >&       aIdsOut,
                    const bool&                   aIsBuildGraph = true ) override
                    {
                        MORIS_ERROR( false, "translate_ids_to_free_point_ids() function has no child implementation" );
                    }

            void translate_ids_to_free_point_ids(
                    const Vector< sint >& aIdsIn,
                    Vector< sint >&       aIdsOut,
                    bool                  aIsBuildGraph = true ) override
                    {
                        MORIS_ERROR( false, "translate_ids_to_free_point_ids() function has no child implementation" );
                    }   

                    // ---------------------------------------------------------------------------------------------------------

                        Matrix< DDSMat >
                map_from_moris_ids_to_indices( const Matrix< DDSMat >& aGlobalIds );

    };
}
