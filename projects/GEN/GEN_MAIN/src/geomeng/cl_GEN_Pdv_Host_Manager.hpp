/*
 * cl_GEN_Pdv_Host_Manager.hpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_

// GEN_MAIN
#include "cl_GEN_Pdv_Host.hpp"
#include "cl_GEN_Field.hpp"
// GEN_CORE
#include "cl_GEN_Dv_Enums.hpp"

// CORE
#include "cl_Matrix.hpp"
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
    // list of pdv hosts - interpolation nodes
    moris::Cell< std::shared_ptr< GEN_Pdv_Host > > mIpHostList;
    // list of pdv hosts - integration nodes
    moris::Cell< std::shared_ptr< GEN_Pdv_Host > > mIgHostList;

    // position in map corresponds to the value of the pdv enum
    Matrix< IndexMat > mIpGlobalPdvTypeMap;
    // position in map corresponds to the value of the pdv enum
    Matrix< IndexMat > mIgGlobalPdvTypeMap;

    // list containing all the used unique pdv types
    moris::Cell< enum GEN_DV > mIpPdvTypeList;
    // list containing all the unique pdv types which are not changing
    moris::Cell< enum GEN_DV > mIpUnchangingPdvTypeList;

    // list containing all the used unique pdv types
    moris::Cell< enum GEN_DV > mIgPdvTypeList;
    // list containing all the unique pdv types which are not changing
    moris::Cell< enum GEN_DV > mIgUnchangingPdvTypeList;

    // total number of dv types on interpolation nodes
    uint mIpNumPDVs = 0;
    // total number of dv types on integration nodes
    uint mIgNumPDVs = 0;

    // global dv ID for interpolation nodes
    uint mIpGlobalID = 0;
    // global dv ID for integration nodes
    uint mIgGlobalID = 0;

//------------------------------------------------------------------------------
public:

//------------------------------------------------------------------------------
    /**
     * constructor
     */
    Pdv_Host_Manager()
    {
        // set size for IP dv type enum to index map
        mIpGlobalPdvTypeMap.set_size( static_cast<size_t>(GEN_DV::END_ENUM), 1, gNoIndex );
        // set size for IG dv type enum to index map
        mIgGlobalPdvTypeMap.set_size( static_cast<size_t>(GEN_DV::ZCOORD) + 1, 1, gNoIndex );
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
    std::shared_ptr< GEN_Pdv_Host > get_ip_pdv_host( moris_index aVertexIndex )
    {
        return mIpHostList( aVertexIndex );
    }
//------------------------------------------------------------------------------
    /**
     * get pdv host for a vertex
     * @param[ out ] aPdvHost a pdv host for the vertex index
     */
    std::shared_ptr< GEN_Pdv_Host > get_ig_pdv_host( moris_index aVertexIndex )
    {
        return mIgHostList( aVertexIndex );
    }
//------------------------------------------------------------------------------
    /**
     * initialize the list of pdv host
     * @param[ in ] aTotalNumVertices number of vertices in IP mesh
     */
    void initialize_ip_hosts( uint aTotalNumVertices )
    {
        // set size for the list of pdv host
        mIpHostList.resize( aTotalNumVertices, nullptr );
    }
//------------------------------------------------------------------------------
    /**
     * initialize the list of pdv host
     * @param[ in ] aTotalNumVertices number of vertices in IP mesh
     */
    void initialize_ig_hosts( uint aTotalNumVertices )
    {
        // set size for the list of pdv host
        mIgHostList.resize( aTotalNumVertices, nullptr );
    }
//------------------------------------------------------------------------------
    /**
     * create pdv host
     * @param[ in ] aNumPdvs     a number of dv
     * @param[ in ] aVertexIndex a vertex index
     */
    void create_ip_pdv_host( uint        aNumPdvs,
                             moris_index aVertexIndex )
    {
        // if pdv host not assigned yet
        if( mIpHostList( aVertexIndex ) == nullptr )
        {
            // create a pdv host
            mIpHostList( aVertexIndex ) = std::make_shared< GEN_Pdv_Host >( aNumPdvs, aVertexIndex );
        }
    }
//------------------------------------------------------------------------------
    /**
     * create pdv host
     * @param[ in ] aNumPdvs     a number of dv
     * @param[ in ] aVertexIndex a vertex index
     */
    void create_ig_pdv_host( uint        aNumPdvs,
                             moris_index aVertexIndex )
    {
        // if pdv host is not assigned yet
        if( mIgHostList( aVertexIndex ) == nullptr )
        {
            // create a pdv host
            mIgHostList( aVertexIndex ) = std::make_shared< GEN_Pdv_Host >( aNumPdvs, aVertexIndex );
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
        this->create_ip_pdv_host( mIpNumPDVs, aVertexIndex );

        // get the pdv host and create the pdv for dv type
        this->get_ip_pdv_host( aVertexIndex )->create_pdv( aPropertyPointer, aPdvType, mIpGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * assign a GEN Field to pdv type by vertex index
     * @param[ in ] aFieldPointer a GEN Field pointer
     * @param[ in ] aPdvType      a list of dv types
     * @param[ in ] aVertexIndex  a vertex index
     */
    void assign_field_to_pdv_type_by_vertex_index( std::shared_ptr< GEN_Field > aFieldPointer,
                                                   enum GEN_DV                  aPdvType,
                                                   moris_index                  aVertexIndex )
    {
        // create a pdv host for vertex index
        this->create_ip_pdv_host( mIpNumPDVs, aVertexIndex );

        // get the pdv host and create the pdv for dv type
        this->get_ip_pdv_host( aVertexIndex )->create_pdv( aFieldPointer, aPdvType, mIpGlobalPdvTypeMap );
    }

//------------------------------------------------------------------------------
    /**
     * FIXME add a phase
     * get pdv by type and vertex index
     * @param[ in ] aVertexIndex     a vertex index
     * @param[ in ] aPdvType         a list of dv types
     */
    std::shared_ptr< GEN_Pdv > get_ip_pdv_by_type_and_index( moris_index aVertexIndex,
                                                             enum GEN_DV aPdvType )
    {
        return this->get_ip_pdv_host( aVertexIndex )->get_pdv_by_type( aPdvType, mIpGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * FIXME add a phase
     * get pdv by type and vertex index
     * @param[ in ] aVertexIndex     a vertex index
     * @param[ in ] aPdvType         a list of dv types
     */
    std::shared_ptr< GEN_Pdv > get_ig_pdv_by_type_and_index( moris_index aVertexIndex,
                                                             enum GEN_DV aPdvType )
    {
        return this->get_ig_pdv_host( aVertexIndex )->get_pdv_by_type( aPdvType, mIgGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * check for active pdv by type and vertex index
     * @param[ in ] aVertexIndex     a vertex index
     * @param[ in ] aPdvType         a list of dv types
     */
    sint check_ip_for_active_types( moris_index aVertexIndex,
                                 enum GEN_DV aPdvType )
    {
        return this->get_ip_pdv_host( aVertexIndex )->is_active_type( aPdvType, mIpGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * check for active pdv by type and vertex index
     * @param[ in ] aVertexIndex     a vertex index
     * @param[ in ] aPdvType         a list of dv types
     */
    sint check_ig_for_active_types( moris_index aVertexIndex,
                                    enum GEN_DV aPdvType )
    {
        return this->get_ig_pdv_host( aVertexIndex )->is_active_type( aPdvType, mIgGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * set dv types
     * @param[ in ] aPdvTypeList list of dv types
     */
    void set_ip_pdv_types( Cell< enum GEN_DV > aPdvTypeList )
    {
        // communicate dv types
        this->communicate_ip_dv_types( aPdvTypeList );

        // get number of dv types
        uint tNumTypes = mIpPdvTypeList.size();

        // loop over dv types
        for( uint i=0; i<tNumTypes; i++)
        {
            // populate the dv type to index map
            mIpGlobalPdvTypeMap( static_cast< sint >( mIpPdvTypeList( i ) ) ) = mIpNumPDVs;

            // update dv type counter
            mIpNumPDVs++;
        }
    }
//------------------------------------------------------------------------------
    /**
     * set dv types
     * @param[ in ] aPdvTypeList list of dv types
     */
    void set_ig_pdv_types( Cell< enum GEN_DV > aPdvTypeList )
    {
        // communicate dv types
        this->communicate_ig_dv_types( aPdvTypeList );

        // get number of dv types
        uint tNumTypes = mIgPdvTypeList.size();

        // loop over dv types
        for( uint i=0; i<tNumTypes; i++)
        {
            // populate the dv type to index map
            mIgGlobalPdvTypeMap( static_cast< sint >( mIgPdvTypeList( i ) ) ) = mIgNumPDVs;

            // update dv type counter
            mIgNumPDVs++;
        }
    }
//------------------------------------------------------------------------------
    /**
     * get pdv type lists
     * @param[ out ] aPdvTypeList list of dv types
     */
    moris::Cell< enum GEN_DV > get_ip_pdv_type_list()
    {
        return mIpPdvTypeList;
    }
//------------------------------------------------------------------------------
    /**
     * get pdv type lists
     * @param[ out ] aPdvTypeList list of dv types
     */
    moris::Cell< enum GEN_DV > get_ig_pdv_type_list()
    {
        return mIgPdvTypeList;
    }
//------------------------------------------------------------------------------
    /**
     * get unchanging pdv type lists
     * @param[ out ] aPdvTypeList list of unchanging dv types
     */
    moris::Cell< enum GEN_DV > get_ip_unchanging_type_list()
    {
        return mIpUnchangingPdvTypeList;
    }
//------------------------------------------------------------------------------
    /**
     * get unchanging pdv type lists
     * @param[ out ] aPdvTypeList list of unchanging dv types
     */
    moris::Cell< enum GEN_DV > get_ig_unchanging_type_list()
    {
        return mIgUnchangingPdvTypeList;
    }
//------------------------------------------------------------------------------
    /**
     * comminucate dv types
     * @param[ in ] aPdvTypeList a local list of dv types
     */
    void communicate_ip_dv_types( moris::Cell< enum GEN_DV > & aPdvTypeList )
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
            mIpPdvTypeList.resize( tNumMaxGlobalDvTypes );
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
                       (mIpPdvTypeList.data()).data(),
                       (tNumLocalDvTypesList.data()).data(),
                       (tDofTypeOffset.data()).data(),
                       MPI_UNSIGNED,
                       0,
                       MPI_COMM_WORLD );

        // temporary variable for mIpPdvTypeList size
        moris::uint tPdvTypeListSize;

        if ( par_rank() == 0 )
        {
            // sort this created list
            std::sort( ( mIpPdvTypeList.data() ).data(), ( mIpPdvTypeList.data() ).data() + mIpPdvTypeList.size() );

            // use std::unique and std::distance to create list containing all used dof types. This list is unique
            auto last = std::unique( ( mIpPdvTypeList.data() ).data(), ( mIpPdvTypeList.data() ).data() + mIpPdvTypeList.size() );
            auto pos  = std::distance( ( mIpPdvTypeList.data() ).data(), last );

            mIpPdvTypeList.resize( pos );

            tPdvTypeListSize = mIpPdvTypeList.size();
        }

        // bcast size of mIpPdvTypeList on processor 0
        MPI_Bcast( & tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

        // resize mIpPdvTypeList on all processors
        mIpPdvTypeList.resize( tPdvTypeListSize );

        // bcast unique mIpPdvTypeList to all processors
        MPI_Bcast( (mIpPdvTypeList.data()).data(), mIpPdvTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    }
//------------------------------------------------------------------------------
    /**
     * comminucate dv types
     * @param[ in ] aPdvTypeList a local list of dv types
     */
    void communicate_ig_dv_types( moris::Cell< enum GEN_DV > & aPdvTypeList )
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
            mIgPdvTypeList.resize( tNumMaxGlobalDvTypes );
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
                       (mIgPdvTypeList.data()).data(),
                       (tNumLocalDvTypesList.data()).data(),
                       (tDofTypeOffset.data()).data(),
                       MPI_UNSIGNED,
                       0,
                       MPI_COMM_WORLD );

        // temporary variable for mIpPdvTypeList size
        moris::uint tPdvTypeListSize;

        if ( par_rank() == 0 )
        {
            // sort this created list
            std::sort( ( mIgPdvTypeList.data() ).data(), ( mIgPdvTypeList.data() ).data() + mIgPdvTypeList.size() );

            // use std::unique and std::distance to create list containing all used dof types. This list is unique
            auto last = std::unique( ( mIgPdvTypeList.data() ).data(), ( mIgPdvTypeList.data() ).data() + mIgPdvTypeList.size() );
            auto pos  = std::distance( ( mIgPdvTypeList.data() ).data(), last );

            mIgPdvTypeList.resize( pos );

            tPdvTypeListSize = mIgPdvTypeList.size();
        }

        // bcast size of mIpPdvTypeList on processor 0
        MPI_Bcast( & tPdvTypeListSize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

        // resize mIpPdvTypeList on all processors
        mIgPdvTypeList.resize( tPdvTypeListSize );

        // bcast unique mIpPdvTypeList to all processors
        MPI_Bcast( (mIgPdvTypeList.data()).data(), mIgPdvTypeList.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD );
    }
//------------------------------------------------------------------------------
    /**
     * update local to global dv type map for interpolation nodes
     */
    void update_ip_local_to_global_dv_type_map(  )
    {
        // get the number of dv types
        uint tNumDvTypes = mIpPdvTypeList.size();

        // get the number of pdv hosts
        uint tNumHosts   = mIpHostList.size();

        // loop over the dv types
        for( uint iType=0; iType<tNumDvTypes; iType++ )
        {
            // loop over the pdv hosts
            for( uint iHost=0; iHost<tNumHosts; iHost++ )
            {
                // if the pdv host was assigned
                if(mIpHostList(iHost) != nullptr)
                {
                    // assign a global id to type
                    mIpHostList(iHost)->assign_id_to_type( mIpGlobalID, mIpPdvTypeList(iType), mIpGlobalPdvTypeMap );

                    // update global id counter
                    mIpGlobalID++;
                }
            }
        }
    }
//------------------------------------------------------------------------------
    /**
     * update local to global dv type map for interpolation nodes
     */
    void update_ig_local_to_global_dv_type_map(  )
    {
        // get the number of dv types
        uint tNumDvTypes = mIgPdvTypeList.size();

        // get the number of pdv hosts
        uint tNumHosts   = mIgHostList.size();

        // loop over the dv types
        for( uint iType=0; iType<tNumDvTypes; iType++ )
        {
            // loop over the pdv hosts
            for( uint iHost=0; iHost<tNumHosts; iHost++ )
            {
                // if the pdv host was assigned
                if(mIgHostList(iHost) != nullptr)
                {
                    // assign a global id to type
                    mIgHostList(iHost)->assign_id_to_type( mIgGlobalID, mIgPdvTypeList(iType), mIgGlobalPdvTypeMap );

                    // update global id counter
                    mIgGlobalID++;
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
    uint get_ip_global_index_for_dv_type( moris_index aVertexIndex,
                                          enum GEN_DV aPdvType )
    {
        return this->get_ip_pdv_host( aVertexIndex )->get_global_index_for_dv_type( aPdvType, mIpGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * get global index for dv type
     * @param[ in ] aVertexIndex a vertex index
     * @param[ in ] aPdvType     a list of dv types
     */
    uint get_ig_global_index_for_dv_type( moris_index aVertexIndex,
                                          enum GEN_DV aPdvType )
    {
        return this->get_ig_pdv_host( aVertexIndex )->get_global_index_for_dv_type( aPdvType, mIgGlobalPdvTypeMap );
    }
//------------------------------------------------------------------------------
    /**
     * get global map for IP nodes
     * @param[ out ] mIpGlobalPdvTypeMap a map from dv type to index
     */
    Matrix< IndexMat > get_ip_global_map(  )
    {
        return mIpGlobalPdvTypeMap;
    }
//------------------------------------------------------------------------------
    /**
     * get global map for IG nodes
     * @param[ out ] mIgGlobalPdvTypeMap a map from dv type to index
     */
    Matrix< IndexMat > get_ig_global_map(  )
    {
        return mIgGlobalPdvTypeMap;
    }
//------------------------------------------------------------------------------
    void mark_ip_pdv_as_unchanging( enum GEN_DV aPdvType )
    {
        mIpUnchangingPdvTypeList.push_back( aPdvType );
    }
//------------------------------------------------------------------------------
    void mark_ig_pdv_as_unchanging( enum GEN_DV aPdvType )
    {
        mIgUnchangingPdvTypeList.push_back( aPdvType );
    }
//------------------------------------------------------------------------------








//--------------- these are temporary for debugging/testing purposes -----------
    /**
     * get_num_hosts
     * @param[ out ] aNumPdvHosts a number of pdv hosts
     */
    uint get_num_ip_hosts()
    {
        return mIpHostList.size();
    }
//------------------------------------------------------------------------------
    uint get_num_ig_hosts()
    {
        return mIgHostList.size();
    }
//------------------------------------------------------------------------------
};
}   // end ge namespace
}       // end moris namepspace



#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_PDV_HOST_MANAGER_HPP_ */
