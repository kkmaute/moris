/*
 * cl_GEN_Design_Variable_Interface.hpp
 *
 *  Created on: Jan 15, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_
#define PROJECTS_GEN_SRC_GEOMENG_CL_GEN_DESIGN_VARIABLE_INTERFACE_HPP_

#include "cl_MSI_Design_Variable_Interface.hpp"

namespace moris
{
namespace ge
{
    class Pdv_Host_Manager;

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
                                     moris::Cell< moris::Matrix< DDSMat > > & aIsActiveDv );
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
                                     moris::Cell< moris::Matrix< DDSMat > > & aIsActiveDv );
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
                                       moris::Matrix< DDRMat >                & aReshapedPdvValues );
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
        Matrix< DDSMat > get_ip_local_global_map();
//------------------------------------------------------------------------------
        /**
         * @brief return local to global dv type map
         */
        Matrix< DDSMat > get_ig_local_global_map();
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
                                                   moris::Cell< moris::Matrix< IdMat > > & aDvIds );
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
                                                   moris::Cell< moris::Matrix< IdMat > > & aDvIds );
//------------------------------------------------------------------------------
        /*
         * @brief returns a cell of all GEN_DV types which are "changing" on the IP mesh
         */
        void get_ip_requested_dv_types( Cell< enum GEN_DV > & aDvTypes );
//------------------------------------------------------------------------------
        /*
         * @brief returns a cell of all GEN_DV types which are "changing" on the IG mesh
         */
        void get_ig_requested_dv_types( Cell< enum GEN_DV > & aDvTypes );
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
