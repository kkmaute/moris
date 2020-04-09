/*
 * cl_GEN_Design_Variable_Interface.hpp
 *
 *  Created on: Jan 15, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_

#include "cl_GEN_Pdv_Host_Manager.hpp"

#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
namespace ge
{
    class GEN_Design_Variable_Interface : MSI::Design_Variable_Interface
    {

    private:

        // pdv host manager pointer
        Pdv_Host_Manager* mPdvHostManager;

        // assembly map for QI
        moris::Cell< moris::Cell< moris_index > > mQIAssemblyMap;

        Matrix< DDSMat > mDummyMat;

    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         * @param[ in ] aPdvHostManager a pdv host manager pointer
         */
        GEN_Design_Variable_Interface( Pdv_Host_Manager* aPdvHostManager )
        : mPdvHostManager( aPdvHostManager )
        {};

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~GEN_Design_Variable_Interface(){};

//------------------------------------------------------------------------------
        /**
         * get dv types for set
         * @param[ in ] aIGMeshSetIndex integration mesh index
         * @param[ in ] aDvTypes        list of groups of dv types to fill
         */
        void get_ip_dv_types_for_set
        ( const moris::moris_index                          aIGMeshSetIndex,
                moris::Cell< moris::Cell< enum GEN_DV > > & aDvTypes )
        {
            // this may change, hold off for now...
            MORIS_ASSERT( false, "GEN_Design_Variable_Interface::get_dv_types_for_set - not implemented." );
        }
//------------------------------------------------------------------------------
        /**
         * get dv types for set
         * @param[ in ] aIGMeshSetIndex integration mesh index
         * @param[ in ] aDvTypes        list of groups of dv types to fill
         */
        void get_ig_dv_types_for_set
        ( const moris::moris_index                          aIGMeshSetIndex,
                moris::Cell< moris::Cell< enum GEN_DV > > & aDvTypes )
        {
            // this may change, hold off for now...
            MORIS_ASSERT( false, "GEN_Design_Variable_Interface::get_dv_types_for_set - not implemented." );
        }
//------------------------------------------------------------------------------
        /**
         * get unique dv types for set
         * @param[ in ] aIGMeshSetIndex integration mesh index
         * @param[ in ] aDvTypes        list dv types to fill
         */
        void get_ip_unique_dv_types_for_set
        ( const moris::moris_index           aIGMeshSetIndex,
                moris::Cell< enum GEN_DV > & aDvTypes )
        {
            MORIS_ASSERT( false, "GEN_Design_Variable_Interface::get_unique_dv_types_for_set - not implemented." );
        }
//------------------------------------------------------------------------------
        /**
         * get unique dv types for set
         * @param[ in ] aIGMeshSetIndex integration mesh index
         * @param[ in ] aDvTypes        list dv types to fill
         */
        void get_ig_unique_dv_types_for_set
        ( const moris::moris_index           aIGMeshSetIndex,
                moris::Cell< enum GEN_DV > & aDvTypes )
        {
            MORIS_ASSERT( false, "GEN_Design_Variable_Interface::get_unique_dv_types_for_set - not implemented." );
        }
//------------------------------------------------------------------------------
        /**
         * get pdv values for requested vertex indices and dv types
         * @param[ in ] aNodeIndices list of node indices
         * @param[ in ] aDvType      list of dv types
         * @param[ in ] aDvValues    list of dv values (DvType)(vertexIndex)
         * @param[ in ] aIsActive    list of active design variables (vertexIndex)(DvType)
         */
        void get_ip_pdv_value( const Matrix< IndexMat >                     & aNodeIndices,
                               const moris::Cell< enum GEN_DV >             & aDvTypes,
                                     moris::Cell< moris::Matrix< DDRMat > > & aDvValues,
                                     moris::Cell< moris::Matrix< DDSMat > > & aIsActiveDv )
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes   = aDvTypes.size();

            // set size for list of dv values
            aDvValues.resize( tNumTypes );

            // set size for list of active flags
            aIsActiveDv.resize( tNumIndices );

            // loop over the node indices
            for( uint iInd =0; iInd < tNumIndices; iInd++ )
            {
                //
                aIsActiveDv( iInd ).resize( tNumTypes, 1 );

                // loop over the requested dv types
                for( uint iType=0; iType < tNumTypes; iType++ )
                {
                    aDvValues( iType ).resize( tNumIndices, 1 );

                    aIsActiveDv( iInd )( iType ) = mPdvHostManager->check_ip_for_active_types( aNodeIndices(iInd), aDvTypes(iType) );

                    if( aIsActiveDv( iInd )( iType ) == 1 )
                    {
                        // TODO: rather than go down the line manager->host->pdv, get the value from a list of pdv values which is stored at the manager level

                        aDvValues( iType )( iInd ) = mPdvHostManager->get_ip_pdv_by_type_and_index( aNodeIndices(iInd), aDvTypes(iType) )->get_val()( 0, 0 );
                    }
                }
            }
        }
//------------------------------------------------------------------------------
        /**
         * get pdv values for requested vertex indices and dv types
         * @param[ in ] aNodeIndices list of node indices
         * @param[ in ] aDvType      list of dv types
         * @param[ in ] aDvValues    list of dv values (DvType)(vertexIndex)
         */
        void get_ip_pdv_value( const Matrix< IndexMat >                     & aNodeIndices,
                               const moris::Cell< enum GEN_DV >             & aDvTypes,
                                     moris::Cell< moris::Matrix< DDRMat > > & aDvValues )
        {
            MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value - overload not implemented on GE side.");
        }
//------------------------------------------------------------------------------
        /**
         * get pdv values for requested vertex indices and dv types
         * @param[ in ] aNodeIndices list of node indices
         * @param[ in ] aDvType      list of dv types
         * @param[ in ] aDvValues    list of dv values (DvType)(vertexIndex)
         * @param[ in ] aIsActive    list of active design variables (vertexIndex)(DvType)
         */
        void get_ig_pdv_value( const Matrix< IndexMat >                     & aNodeIndices,
                               const moris::Cell< enum GEN_DV >             & aDvTypes,
                                     moris::Cell< moris::Matrix< DDRMat > > & aDvValues,
                                     moris::Cell< moris::Matrix< DDSMat > > & aIsActiveDv )
        {
            // get the number of node indices requested
            uint tNumIndices = aNodeIndices.length();

            // get the number of dv types requested
            uint tNumTypes   = aDvTypes.size();

            // set size for list of dv values
            aDvValues.resize( tNumTypes );

            // set size for list of active flags
            aIsActiveDv.resize( tNumIndices );

            // loop over the node indices
            for( uint iInd =0; iInd < tNumIndices; iInd++ )
            {
                //
                aIsActiveDv( iInd ).resize( tNumTypes, 1 );

                // loop over the requested dv types
                for( uint iType=0; iType < tNumTypes; iType++ )
                {
                    aDvValues( iType ).resize( tNumIndices, 1 );

                    aIsActiveDv( iInd )( iType ) = mPdvHostManager->check_ig_for_active_types( aNodeIndices(iInd), aDvTypes(iType) );

                    if( aIsActiveDv( iInd )( iType ) == 1 )
                    {
                        // TODO: rather than go down the line manager->host->pdv, get the value from a list of pdv values which is stored at the manager level
                        aDvValues( iType )( iInd ) = mPdvHostManager->get_ig_pdv_by_type_and_index( aNodeIndices(iInd), aDvTypes(iType) )->get_val()( 0, 0 );
                    }
                }
            }
        }
//------------------------------------------------------------------------------
        /**
         * get pdv values for requested vertex indices and dv types
         * @param[ in ] aNodeIndices list of node indices
         * @param[ in ] aDvType      list of dv types
         * @param[ in ] aDvValues    list of dv values (DvType)(vertexIndex)
         */
        void get_ig_pdv_value( const Matrix< IndexMat >                     & aNodeIndices,
                               const moris::Cell< enum GEN_DV >             & aDvTypes,
                                     moris::Cell< moris::Matrix< DDRMat > > & aDvValues )
        {
            MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value - overload not implemented on GE side.");
        }

//------------------------------------------------------------------------------
        /**
         * reshape pdv values
         * i.e. reshape a cell of matrix to a matrix
         * @param[ in ] aPdvValues         a cell of matrices with pdv values
         * @param[ in ] aReshapedPdvValues a matrix of pdv values
         */
        void reshape_pdv_values( const moris::Cell< moris::Matrix< DDRMat > > & aPdvValues,
                                       moris::Matrix< DDRMat >                & aReshapedPdvValues )
        {
            MORIS_ASSERT( aPdvValues.size() != 0,
                          "GEN_Design_Variable_Interface::reshape_pdv_value - pdv value vector is empty.");

            // get the number of rows and columns
            uint tRows = aPdvValues( 0 ).numel();
            uint tCols = aPdvValues.size();

            // set size for the reshaped matrix
            aReshapedPdvValues.set_size( tRows, tCols );

            for( uint iCol = 0; iCol < tCols; iCol++ )
            {
                aReshapedPdvValues( { 0, tRows - 1 }, { iCol, iCol } )
                = aPdvValues( iCol ).matrix_data();
            }
        }

//------------------------------------------------------------------------------
        /**
         * return local to global dv map.
         * ( this is a collection of all local parallel consistent Ids )
         */
        moris::Matrix< DDSMat > get_local_global_map()
        {
            MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value - overload not implemented on GE side.");
            return mDummyMat;
        };

//------------------------------------------------------------------------------
        /**
         * return owned local to global dv map.
         * ( this is a collection of all owned local parallel consistent Ids )
         */
        moris::Matrix< DDSMat > get_owned_local_global_map()
        {
            MORIS_ASSERT(false, "Design_Variable_Interface::get_pdv_value - overload not implemented on GE side.");
            return mDummyMat;
        };

//------------------------------------------------------------------------------
        /**
         * @brief return local to global dv type map
         */
        Matrix< DDSMat > get_ip_local_global_map()
        {
            return mPdvHostManager->get_ip_global_map();
        }
//------------------------------------------------------------------------------
        /**
         * @brief return local to global dv type map
         */
        Matrix< DDSMat > get_ig_local_global_map()
        {
            return mPdvHostManager->get_ig_global_map();
        }
//------------------------------------------------------------------------------
        /**
         * @brief return local to global DV type map
         * @param[ in ] aVertexIndex   List of vertex indices
         * @param[ in ] aDvType        List of Dv types
         * @param[ in ] aDvIds         List of Dv Ids
         *
         */
        void get_ip_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index >     & aNodeIndices,
                                             const moris::Cell< enum GEN_DV >            & aDvTypes,
                                                   moris::Cell< moris::Matrix< IdMat > > & aDvIds )
        {
            /*
             * - each cell is a row vector of global IDs per each type
             * - return the global ids of the dv type on a specified vertex
             */

            uint tNumIndices = aNodeIndices.size();
            uint tNumTypes   = aDvTypes.size();

            moris::Cell< uint > tCounter( tNumTypes, 0 );

            aDvIds.resize( tNumTypes );
            for( uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds( Ik ).set_size( tNumIndices, 1);
            }

            for( uint iType=0; iType<tNumTypes; iType++ )
            {
                for( uint iInd=0; iInd<tNumIndices; iInd++ )
                {

                    bool tDvTypeExists = mPdvHostManager->check_ip_for_active_types( aNodeIndices(iInd), aDvTypes(iType) );   // flag for if the DV type exists on the current host

                    if( tDvTypeExists )
                    {
                        aDvIds( iType )( iInd ) = mPdvHostManager->get_ip_global_index_for_dv_type( aNodeIndices(iInd), aDvTypes(iType) );
                        tCounter( iType )++;
                    }
                }
            }

            for( uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds( Ik ).resize( tCounter( Ik ), 1);
            }
        }
//------------------------------------------------------------------------------
        /**
         * @brief return local to global DV type map
         * @param[ in ] aVertexIndex   List of vertex indices
         * @param[ in ] aDvType        List of Dv types
         * @param[ in ] aDvIds         List of Dv Ids
         *
         */
        void get_ig_dv_ids_for_type_and_ind( const moris::Cell< moris::moris_index >     & aNodeIndices,
                                             const moris::Cell< enum GEN_DV >            & aDvTypes,
                                                   moris::Cell< moris::Matrix< IdMat > > & aDvIds )
        {
            /*
             * - each cell is a row vector of global IDs per each type
             * - return the global ids of the dv type on a specified vertex
             */

            uint tNumIndices = aNodeIndices.size();
            uint tNumTypes   = aDvTypes.size();

            moris::Cell< uint > tCounter( tNumTypes, 0 );

            aDvIds.resize( tNumTypes );
            for( uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds( Ik ).set_size( tNumIndices, 1);
            }

            for( uint iType=0; iType<tNumTypes; iType++ )
            {
                for( uint iInd=0; iInd<tNumIndices; iInd++ )
                {

                    bool tDvTypeExists = mPdvHostManager->check_ig_for_active_types( aNodeIndices(iInd), aDvTypes(iType) );   // flag for if the DV type exists on the current host

                    if( tDvTypeExists )
                    {
                        aDvIds( iType )( iInd ) = mPdvHostManager->get_ig_global_index_for_dv_type( aNodeIndices(iInd), aDvTypes(iType) );
                        tCounter( iType )++;
                    }
                }
            }

            for( uint Ik = 0; Ik <tNumTypes; Ik++)
            {
                aDvIds( Ik ).resize( tCounter( Ik ), 1);
            }
        }
//------------------------------------------------------------------------------
        /*
         * @brief returns a cell of all GEN_DV types which are "changing" on the IP mesh
         */
        void get_ip_requested_dv_types( Cell< enum GEN_DV > & aDvTypes )
        {
            // total number of unique DV types
            uint tNumDvTypes = mPdvHostManager->get_ip_pdv_type_list().size();

            // list of all DV types which are changing
            aDvTypes.reserve(tNumDvTypes);

            // total number of unchanging DV types
            uint tNumUnchangingDvTypes = mPdvHostManager->get_ip_unchanging_type_list().size();

            for( uint i=0; i<tNumDvTypes; i++ )
            {
                bool tIsChanging;
                for( uint j=0; j<tNumUnchangingDvTypes; j++ )
                {
                    if( mPdvHostManager->get_ip_pdv_type_list()(i) == mPdvHostManager->get_ip_unchanging_type_list()(j) )
                    {
                        // this type is not changing and therefore is not added to the list
                        tIsChanging = false;
                        break;
                    }
                    else
                    {
                        // this type is changing and so we add it to the list
                        tIsChanging = true;
                    }
                }
                // use flag to make assignment
                if(tIsChanging)
                {
                    aDvTypes.push_back( mPdvHostManager->get_ip_pdv_type_list()(i) );
                }
            }

            // shrink list to appropriate size
            aDvTypes.shrink_to_fit();

            // get rid of redundant entries
            unique(aDvTypes);
        }
//------------------------------------------------------------------------------
        /*
         * @brief returns a cell of all GEN_DV types which are "changing" on the IG mesh
         */
        void get_ig_requested_dv_types( Cell< enum GEN_DV > & aDvTypes )
        {
            // total number of unique DV types
            uint tNumDvTypes = mPdvHostManager->get_ig_pdv_type_list().size();

            // list of all DV types which are changing
            aDvTypes.reserve(tNumDvTypes);

            // total number of unchanging DV types
            uint tNumUnchangingDvTypes = mPdvHostManager->get_ig_unchanging_type_list().size();

            for( uint i=0; i<tNumDvTypes; i++ )
            {
                bool tIsChanging;
                for( uint j=0; j<tNumUnchangingDvTypes; j++ )
                {
                    if( mPdvHostManager->get_ig_pdv_type_list()(i) == mPdvHostManager->get_ig_unchanging_type_list()(j) )
                    {
                        // this type is not changing and therefore is not added to the list
                        tIsChanging = false;
                        break;
                    }
                    else
                    {
                        // this type is changing and so we add it to the list
                        tIsChanging = true;
                    }
                }
                // use flag to make assignment
                if(tIsChanging)
                {
                    aDvTypes.push_back( mPdvHostManager->get_ig_pdv_type_list()(i) );
                }
            }

            // shrink list to appropriate size
            aDvTypes.shrink_to_fit();

            // get rid of redundant entries
            unique(aDvTypes);
        }
//------------------------------------------------------------------------------
        /**
         * get QI assembly map
         */
        moris::Cell< moris::Cell< moris_index > > & get_QI_assembly_map()
        {
            return mQIAssemblyMap;
        }
    };

}   // end ge namespace
}   // end moris namespace

#endif /* PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_ */
