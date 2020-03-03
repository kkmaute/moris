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
    // list of pdv host
    moris::Cell< std::shared_ptr< GEN_Pdv_Host > > mHostList;

    // position in map corresponds to the value of the pdv enum
    Matrix< IndexMat > mGlobalPdvTypeMap;

    // list containing all the used unique pdv types
    moris::Cell< enum GEN_DV > mPdvTypeList;

    // total number of dv types
    uint mNumPDVs = 0;

    // global dv ID
    uint mGlobalID = 0;

//------------------------------------------------------------------------------
public:

//------------------------------------------------------------------------------
    /**
     * constructor
     */
    Pdv_Host_Manager()
    {
        // set size for dv type enum to index map
        mGlobalPdvTypeMap.set_size( static_cast<size_t>(GEN_DV::END_ENUM), 1, gNoIndex );
    };

//------------------------------------------------------------------------------
    /**
     * trivial destructor
     */
    ~Pdv_Host_Manager(){};

//------------------------------------------------------------------------------
    /**
     * get pdv host for a vertex
     * @param[ out ] aPdvHost a pdv host for the vertex index
     */
    std::shared_ptr< GEN_Pdv_Host > get_pdv_host( moris_index aVertexIndex )
    {
        return mHostList( aVertexIndex );
    }

//------------------------------------------------------------------------------
    /**
     * initialize the list of pdv host
     * @param[ in ] aTotalNumVertices number of vertices in IP mesh
     */
    void initalize_hosts( uint aTotalNumVertices )
    {
        // set size for the list of pdv host
        mHostList.resize( aTotalNumVertices, nullptr );
    }

//    //------------------------------------------------------------------------------
//    void initialize_new_hosts( uint aNumNewNodes )    // assuming the indices are consecutive
//    {
//        for( uint i=0; i<aNumNewNodes; i++ )
//        {
//            mHostList.push_back(nullptr);
//        }
//    }

//------------------------------------------------------------------------------
    /**
     * create pdv host
     * @param[ in ] aNumPdvs     a number of dv
     * @param[ in ] aVertexIndex a vertex index
     */
    void create_pdv_host( uint        aNumPdvs,
                          moris_index aVertexIndex )
    {
        // if pdv host not assigned yet
        if( mHostList( aVertexIndex ) == nullptr )
        {
            // create a pdv host
            mHostList( aVertexIndex ) = std::make_shared< GEN_Pdv_Host >( aNumPdvs, aVertexIndex );
        }
    }

//------------------------------------------------------------------------------
    /**
     * assign a GEN property to pdv type by vertex index
     * @param[ in ] aPropertyPointer a GEN property pointer
     * @param[ in ] aPdvType         a list of dv types
     * @param[ in ] aVertexIndex     a vertex index
     */
    void assign_property_to_pdv_type_by_vertex_index( std::shared_ptr< GEN_Property > aPropertyPointer,
                                                      enum GEN_DV                     aPdvType,
                                                      moris_index                     aVertexIndex )
    {
        // create a pdv host for vertex index
        this->create_pdv_host( mNumPDVs, aVertexIndex );

        // get the pdv host and create the pdv for dv type
        this->get_pdv_host( aVertexIndex )->create_pdv( aPropertyPointer, aPdvType, mGlobalPdvTypeMap );
    }

//------------------------------------------------------------------------------
    /**
     * FIXME add a phase
     * get pdv by type and vertex index
     * @param[ in ] aVertexIndex     a vertex index
     * @param[ in ] aPdvType         a list of dv types
     */
    std::shared_ptr< GEN_Pdv > get_pdv_by_type_and_index( moris_index aVertexIndex,
                                                          enum GEN_DV aPdvType )
    {
        return this->get_pdv_host( aVertexIndex )->get_pdv_by_type( aPdvType, mGlobalPdvTypeMap );
    }

//------------------------------------------------------------------------------
    /**
     * check for active pdv by type and vertex index
     * @param[ in ] aVertexIndex     a vertex index
     * @param[ in ] aPdvType         a list of dv types
     */
    sint check_for_active_types( moris_index aVertexIndex,
                                 enum GEN_DV aPdvType )
    {
        return this->get_pdv_host( aVertexIndex )->is_active_type( aPdvType, mGlobalPdvTypeMap );
    }

//------------------------------------------------------------------------------
    /**
     * set dv types
     * @param[ in ] aPdvTypeList list of dv types
     */
    void set_pdv_types( Cell< enum GEN_DV > aPdvTypeList )
    {
        // communicate dv types
        this->communicate_dv_types( aPdvTypeList );

        // get number of dv types
        uint tNumTypes = mPdvTypeList.size();

        // loop over dv types
        for( uint i=0; i<tNumTypes; i++)
        {
            // populate the dv type to index map
            mGlobalPdvTypeMap( static_cast< sint >( mPdvTypeList( i ) ) ) = mNumPDVs;

            // update dv type counter
            mNumPDVs++;
        }
    }

//------------------------------------------------------------------------------
    /**
     * get pdv type lists
     * @param[ out ] aPdvTypeList list of dv types
     */
    moris::Cell< enum GEN_DV > get_pdv_type_list()
    {
        return mPdvTypeList;
    }

//------------------------------------------------------------------------------
    /**
     * comminucate dv types
     * @param[ in ] aPdvTypeList a local list of dv types
     */
    void communicate_dv_types( moris::Cell< enum GEN_DV > & aPdvTypeList )
    {
        // get processor size
        int tSize = par_size();

        // get number of local dv types
        moris::sint tNumLocalDvTypes = aPdvTypeList.size();

        // variable for maximal possible global dv types ???
        moris::sint tNumMaxGlobalDvTypes;

        // get number of global dv types
        sum_all( tNumLocalDvTypes, tNumMaxGlobalDvTypes );

        if ( par_rank() == 0 )
        {
            // set size of dv type list = number of global types
            mPdvTypeList.resize( tNumMaxGlobalDvTypes );
        }

        // create list containing the number of local dv types
        moris::Cell < moris::sint > tNumLocalDvTypesList ( tSize );

        // insert number of local dv types into list containing the number of local dv types
        MPI_Allgather( &tNumLocalDvTypes, 1, MPI_UNSIGNED, (tNumLocalDvTypesList.data()).data(), 1, MPI_UNSIGNED,  MPI_COMM_WORLD );

        // create list containing the offsets of the local dv types in relation to processor 0
        moris::Cell< moris::sint > tDofTypeOffset( tSize, 0 );

        // fill the list with the corresponding offsets
        for ( int Ip = 1; Ip < tSize; ++Ip )
        {
            tDofTypeOffset( Ip ) = tDofTypeOffset( Ip-1 ) + tNumLocalDvTypesList( Ip-1 );
        }

        // assemble list containing all used dv types. Dv types are not unique
        MPI_Gatherv( ((aPdvTypeList.data()).data()),
                       tNumLocalDvTypes,
                       MPI_UNSIGNED,
                       (mPdvTypeList.data()).data(),
                       (tNumLocalDvTypesList.data()).data(),
                       (tDofTypeOffset.data()).data(),
                       MPI_UNSIGNED,
                       0,
                       MPI_COMM_WORLD );

        // temporary variable for mPdvTypeList size
        moris::uint tPdvTypeListSize;

        if ( par_rank() == 0 )
        {
            // sort this created list
            std::sort( ( mPdvTypeList.data() ).data(), ( mPdvTypeList.data() ).data() + mPdvTypeList.size() );

            // use std::unique and std::distance to create list containing all used dof types. This list is unique
            auto last = std::unique( ( mPdvTypeList.data() ).data(), ( mPdvTypeList.data() ).data() + mPdvTypeList.size() );
            auto pos  = std::distance( ( mPdvTypeList.data() ).data(), last );

            mPdvTypeList.resize( pos );

            tPdvTypeListSize = mPdvTypeList.size();
        }

        // bcast size of mPdvTypeList on processor 0
        MPI_Bcast( & tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

        // resize mPdvTypeList on all processors
        mPdvTypeList.resize( tPdvTypeListSize );

        // bcast unique mPdvTypeList to all processors
        MPI_Bcast( (mPdvTypeList.data()).data(), mPdvTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    }

//------------------------------------------------------------------------------
    /**
     * update local to global dv type map
     */
    void update_local_to_global_dv_type_map()
    {
        // get the number of dv types
        uint tNumDvTypes = mPdvTypeList.size();

        // get the number of pdv hosts
        uint tNumHosts   = mHostList.size();

        // loop over the dv types
        for( uint iType=0; iType<tNumDvTypes; iType++ )
        {
            // loop over the pdv hosts
            for( uint iHost=0; iHost<tNumHosts; iHost++ )
            {
                // if the pdv host was assigned
                if(mHostList(iHost) != nullptr)
                {
                    // assign a global id to type
                    mHostList(iHost)->assign_id_to_type( mGlobalID, mPdvTypeList(iType), mGlobalPdvTypeMap );

                    // update global id counter
                    mGlobalID++;
                }
            }
        }
    }

//------------------------------------------------------------------------------
    /**
     * get global index for dv type
     * @param[ in ] aVertexIndex a vertex index
     * @param[ in ] aPdvType     a list of dv types
     */
    uint get_global_index_for_dv_type( moris_index aVertexIndex,
                                       enum GEN_DV aPdvType )
    {
        return this->get_pdv_host( aVertexIndex )->get_global_index_for_dv_type( aPdvType, mGlobalPdvTypeMap );
    }

//------------------------------------------------------------------------------
    /**
     * get global map
     * @param[ out ] mGlobalPdvTypeMap a map from dv type to index
     */
    Matrix< IndexMat > get_global_map(  )
    {
        return mGlobalPdvTypeMap;
    }

//------------------------------------------------------------------------------
    // these are temporary for debugging/testing purposes
    /**
     * get_num_hosts
     * @param[ out ] aNumPdvHosts a number of pdv hosts
     */
    uint get_num_hosts()
    {
        return mHostList.size();
    }

    //------------------------------------------------------------------------------
};
}   // end ge namespace
}       // end moris namepspace



#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
