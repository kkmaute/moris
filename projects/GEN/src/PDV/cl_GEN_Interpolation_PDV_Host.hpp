/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Interpolation_PDV_Host.hpp
 *
 */

#pragma once

#include "GEN_Data_Types.hpp"
#include "cl_GEN_PDV.hpp"

namespace moris::ge
{
    class Property;
    class PDV_Host_Manager;

    class Interpolation_PDV_Host
    {
      private:
        PDV_Host_Manager* mPDVHostManager = nullptr;

        // Identifies the host node
        moris_index mNodeIndex;
        moris_id    mNodeId;
        moris_index mNodeOwner;

        Matrix< DDRMat > mCoordinates;

        // Information about the contained PDVs
        Cell< std::shared_ptr< PDV > > mPDVs;

      public:
        /**
         * Constructor
         *
         * @param aPDVTypes PDV types for this host
         * @param aStartingGlobalIndex Global index to start assigning new PDV types
         */
        Interpolation_PDV_Host(
                PDV_Host_Manager*       aPDVHostManager,
                const moris_index&      aNodeIndex,
                const moris_id&         aNodeId,
                const moris_index&      aNodeOwner,
                const Matrix< DDRMat >& aCoordinates );

        /**
         * destructor
         */
        ~Interpolation_PDV_Host();

        /**
         * Gets the number of PDVs on this PDV host.
         *
         * @return Number of PDVs
         */
        uint get_num_pdvs();

        /**
         * Gets the Id of this PDV host,
         *
         * @return PDV host id
         */
        moris_id
        get_id()
        {
            return mNodeId;
        };

        /**
         * Gets rank of owning processor,
         *
         * @return rank of owning processor
         */
        moris_index
        get_pdv_owning_processor()
        {
            return mNodeOwner;
        };

        /**
         * Sets PDV ID for given type.
         *
         * @param[in] aPDVType PDV type
         * @param[in] aId      PDV Id
         */
        void set_pdv_id( PDV_Type aPDVType, const moris_id aId );

        /**
         * Gets PDV ID by type.
         *
         * @param aPDVType PDV type
         * @return PDV ID
         */
        moris_id get_pdv_id( PDV_Type aPDVType );

        /**
         * Gets PDV ID by index.
         *
         * @param aPDVIndex PDV index
         * @return PDV ID
         */
        moris_id get_pdv_id( uint aPDVIndex );

        /**
         * Create PDV with real value.
         *
         * @param aPDVType PDV type
         * @param aPDVVal PDV value
         */
        void create_pdv( PDV_Type aPDVType, moris::real aPDVVal );

        /**
         * Create PDV with GEN property.
         *
         * @param aPDVType PDV type
         * @param aPropertyPointer Pointer to a GEN property
         */
        void create_pdv( PDV_Type aPDVType, std::shared_ptr< Property > aPropertyPointer );

        /**
         * Check if PDV type is active on this host.
         *
         * @param aPDVType PDV type
         * @return if PDV type is active
         */
        bool is_active_type( PDV_Type aPDVType );

        /**
         * Get global index for pdv by type.
         *
         * @param aPDVType PDV type
         * @return Global index
         */
        void set_global_index_for_pdv_type( PDV_Type aPDVType, moris_id aId );

        /**
         * Get global index for pdv by type.
         *
         * @param aPDVType PDV type
         * @return Global index
         */
        uint get_global_index_for_pdv_type( PDV_Type aPDVType );

        /**
         * Get all of the global PDV indices on this host.
         *
         * @param aGlobalPDVIndices matrix of indices to be returned
         */
        Matrix< DDSMat > get_all_global_indices();

        /**
         * Get the value of a PDV by type.
         *
         * @param aPDVType PDV type
         * @return Value on this PDV
         */
        real get_pdv_value( PDV_Type aPDVType );

        /**
         * Gets whether or not this PDV exists on this host.
         *
         * @param aPDVType PDV type
         * @return if the type exists
         */
        bool get_pdv_exists( PDV_Type aPDVType );

        /**
         * Gets all of the sensitivity vectors on each PDV.
         *
         * @param aPDVIndex PDV index
         * @return Sensitivities for this PDV index
         */
        Matrix< DDRMat > get_sensitivities( uint aPDVIndex );

        /**
         * Gets the IDs of ADVs which the given PDV depends on.
         *
         * @param aPDVIndex PDV index
         * @return ADV IDs
         */
        Matrix< DDSMat > get_determining_adv_ids( uint aPDVIndex );

        /**
         * Gets coordinates stored with host.
         * @return coordinates
         */
        const Matrix< DDRMat >&
        get_coords()
        {
            return mCoordinates;
        }

        /**
         * Print current state of PDV host
         */
        void print_state();
    };
}
