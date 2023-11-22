/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Adof.hpp
 *
 */

#ifndef SRC_FEM_CL_ADOF_HPP_
#define SRC_FEM_CL_ADOF_HPP_

#include "moris_typedefs.hpp"

namespace moris
{
    namespace MSI
    {
        class Adof
        {
          private:
            moris_id    mAdofId          = gNoID;
            moris_id    mAdofExternalId  = gNoID;       // FIXME delete
            moris_index mAdofExternalInd = gNoIndex;    // external Ind. only for HMR use

            uint mOwningProcessor;

            sint mAdofTypeTimeIdentifier = -1;    // Stores a unique integer for every dof type and time combination.
            // only used for multigrid

          public:
            Adof(){};

            ~Adof(){};

            void
            set_adof_id( const moris_id aAdofId )
            {
                mAdofId = aAdofId;
            };

            void
            set_adof_external_id( const moris_id aAdofExtId )    // FIXME delte
            {
                mAdofExternalId = aAdofExtId;
            };

            void
            set_adof_external_ind( const moris_index aAdofExtInd )
            {
                mAdofExternalInd = aAdofExtInd;
            };

            void
            set_adof_owning_processor( const sint aOwningProcessor )
            {
                mOwningProcessor = aOwningProcessor;
            };

            moris_id
            get_adof_id()
            {
                return mAdofId;
            };

            moris_id
            get_adof_external_id()
            {
                return mAdofExternalId;
            };

            moris_index
            get_adof_external_ind()
            {
                return mAdofExternalInd;
            };

            moris_id
            get_adof_owning_processor()
            {
                return mOwningProcessor;
            };

            sint
            get_adof_type_time_identifier()
            {
                return mAdofTypeTimeIdentifier;
            };

            void
            set_adof_type_time_identifier( const sint aAdofTypeTimeIdentifier )
            {
                mAdofTypeTimeIdentifier = aAdofTypeTimeIdentifier;
            };

        };  // class ADof
    }   // namespace MSI
}    // namespace moris

#endif /* SRC_FEM_CL_ADOF_HPP_ */
