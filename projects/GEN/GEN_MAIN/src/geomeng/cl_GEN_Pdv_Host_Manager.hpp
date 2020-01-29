/*
 * cl_GEN_Pdv_Host_Manager.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_

// GEN
#include "cl_GEN_Pdv_Host.hpp"
#include "cl_Matrix.hpp"

#include "cl_GEN_Dv_Enums.hpp"

// CORE
#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
namespace ge
{
class Pdv_Host_Manager
{
private:
    moris::Cell< std::shared_ptr< GEN_Pdv_Host > > mHostList;

    Matrix< IndexMat > mGlobalPdvTypeMap;      // position in map corresponds to the value of the pdv enum

    moris::Cell< enum GEN_DV > mPdvTypeList;   // list containing all the used unique pdv types

    uint mNumPDVs = 0;                         // total number of dv types

    uint mGlobalID = 0;                        // global dv ID

    //------------------------------------------------------------------------------
public:
    Pdv_Host_Manager( )
    {
        mGlobalPdvTypeMap.set_size( static_cast<size_t>(GEN_DV::END_ENUM), 1, gNoIndex );
    };

    ~Pdv_Host_Manager(){};
    //------------------------------------------------------------------------------



    std::shared_ptr< GEN_Pdv_Host > get_pdv_host( moris_index aVertexIndex )
    {
        return mHostList( aVertexIndex );
    }
    //------------------------------------------------------------------------------
    void initalize_hosts( uint aTotalNumVertices )
    {
        mHostList.resize( aTotalNumVertices, nullptr );
    }
    //------------------------------------------------------------------------------
    void create_pdv_host( uint        aNumPdvs,
                          moris_index aVertexIndex )
    {
        mHostList( aVertexIndex ) = std::make_shared< GEN_Pdv_Host >( aNumPdvs, aVertexIndex );
    }
    //------------------------------------------------------------------------------
    void assign_property_to_pdv_type_by_vertex_index( std::shared_ptr< GEN_Property > aPropertyPointer,
                                                      enum GEN_DV                     aPdvType,
                                                      moris_index                     aVertexIndex )
    {
        this->create_pdv_host( mNumPDVs, aVertexIndex );
        this->get_pdv_host( aVertexIndex )->create_pdv( aPropertyPointer, aPdvType, mGlobalPdvTypeMap );
    }
    //------------------------------------------------------------------------------
    std::shared_ptr< GEN_Pdv > get_pdv_by_type_and_index( moris_index aVertexIndex,
                                                          enum GEN_DV aPdvType )
    {
        return this->get_pdv_host( aVertexIndex )->get_pdv_by_type( aPdvType, mGlobalPdvTypeMap );
    }
    //------------------------------------------------------------------------------
    sint check_for_active_types( moris_index aVertexIndex,
                                 enum GEN_DV aPdvType )
    {
        return this->get_pdv_host( aVertexIndex )->is_active_type( aPdvType, mGlobalPdvTypeMap );
    }
    //------------------------------------------------------------------------------
    void set_pdv_types( Cell< enum GEN_DV > aPdvTypeList )
    {
        uint tNumTypes = aPdvTypeList.size();
        for( uint i=0; i<tNumTypes; i++)
        {
            mGlobalPdvTypeMap( static_cast<sint>(aPdvTypeList(i)) ) = mNumPDVs;
            mNumPDVs++;
        }

        this->communicate_dof_types( aPdvTypeList );
    }
    //------------------------------------------------------------------------------
    moris::Cell< enum GEN_DV > get_pdv_type_list(  )
    {
        return mPdvTypeList;
    }
    //------------------------------------------------------------------------------
    void communicate_dof_types( moris::Cell< enum GEN_DV > & aPdvTypeList )
    {
        // Get processor size
        int tSize = par_size();

        // Get number of local dof types
        moris::sint tNumLocalDofTypes = aPdvTypeList.size();

        // Variable for maximal possible global dof types
        moris::sint tNumMaxGlobalDofTypes;

        // Get number of global dof types
        sum_all( tNumLocalDofTypes, tNumMaxGlobalDofTypes );

        if ( par_rank() == 0 )
        {
            // Set size of of pdof type list = number of global types
            mPdvTypeList.resize( tNumMaxGlobalDofTypes );
        }

        // Create list containing the number of local dof types
        moris::Cell < moris::sint > tNumLocalDofTypesList ( tSize );

        // Insert number of local dof types into list containing the number of local dof types
        MPI_Allgather( &tNumLocalDofTypes, 1, MPI_UNSIGNED, (tNumLocalDofTypesList.data()).data(), 1, MPI_UNSIGNED,  MPI_COMM_WORLD );

        // Create list containing the offsets of the local dof types in relation to processor 0
        moris::Cell< moris::sint > tDofTypeOffset( tSize, 0 );

        // Fill the list with the corresponding offsets
        for ( int Ip = 1; Ip < tSize; ++Ip )
        {
            tDofTypeOffset( Ip ) = tDofTypeOffset( Ip-1 ) + tNumLocalDofTypesList( Ip-1 );
        }

        // Assemble list containing all used dof types. Dof types are not unique
        MPI_Gatherv( ((aPdvTypeList.data()).data()),
                tNumLocalDofTypes,
                MPI_UNSIGNED,
                (mPdvTypeList.data()).data(),
                (tNumLocalDofTypesList.data()).data(),
                (tDofTypeOffset.data()).data(),
                MPI_UNSIGNED,
                0,
                MPI_COMM_WORLD );

        // Temporary variable for mPdvTypeList size
        moris::uint tPdofTypeListSize;

        if ( par_rank() == 0 )
        {
            // Sort this created list
            std::sort( ( mPdvTypeList.data() ).data(), ( mPdvTypeList.data() ).data() + mPdvTypeList.size() );

            // use std::unique and std::distance to create list containing all used dof types. This list is unique
            auto last = std::unique( ( mPdvTypeList.data() ).data(), ( mPdvTypeList.data() ).data() + mPdvTypeList.size() );
            auto pos  = std::distance( ( mPdvTypeList.data() ).data(), last );

            mPdvTypeList.resize( pos );

            tPdofTypeListSize = mPdvTypeList.size();
        }

        // Bcast size of mPdvTypeList on processor 0
        MPI_Bcast( & tPdofTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

        // Resize mPdvTypeList on all processors
        mPdvTypeList.resize( tPdofTypeListSize );

        // Bcast unique mPdvTypeList to all processors
        MPI_Bcast( (mPdvTypeList.data()).data(), mPdvTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    }
    //------------------------------------------------------------------------------
    void update_local_to_global_dv_type_map(  )
    {
        uint tNumDvTypes = mPdvTypeList.size();
        uint tNumHosts   = mHostList.size();

        for( uint iType=0; iType<tNumDvTypes; iType++ )
        {
            for( uint iHost=0; iHost<tNumHosts; iHost++ )
            {
                if(mHostList(iHost) != nullptr)
                {
                    mHostList(iHost)->assign_id_to_type( mGlobalID, mPdvTypeList(iType), mGlobalPdvTypeMap );
                    mGlobalID++;
                }
            }
        }
    }
    //------------------------------------------------------------------------------
    uint get_global_index_for_dv_type( moris_index aVertexIndex, enum GEN_DV aPdvType )
    {
        return this->get_pdv_host( aVertexIndex )->get_global_index_for_dv_type( aPdvType, mGlobalPdvTypeMap );
    }
    //------------------------------------------------------------------------------
    Matrix< IndexMat > get_global_map(  )
    {
        return mGlobalPdvTypeMap;
    }











};
}   // end ge namespace
}       // end ,moris namepspace



#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
