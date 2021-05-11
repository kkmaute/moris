/*
 * cl_Adof.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_ADOF_HPP_
#define SRC_FEM_CL_ADOF_HPP_

#include "typedefs.hpp"

namespace moris
{
    namespace MSI
    {
        class Adof
        {
        private:
            moris::moris_id mAdofId = gNoID;
            moris::moris_id mAdofExternalId = gNoID; //FIXME delete
            moris::moris_index mAdofExternalInd = gNoIndex; //external Ind. only for HMR use

            moris::uint mOwningProcessor;

            moris::sint mAdofTypeTimeIdentifier = -1;           // Stores a unique integer for every dof type and time combination.
            // only used for multigrid

        public:
            Adof()
            {};

            ~Adof()
            {};

            void set_adof_id( const moris::moris_id aAdofId )
            {
                mAdofId = aAdofId;
            };

            void set_adof_external_id( const moris::moris_id aAdofExtId )   //FIXME delte
            {
                mAdofExternalId = aAdofExtId;
            };

            void set_adof_external_ind( const moris::moris_index aAdofExtInd )
            {
                mAdofExternalInd = aAdofExtInd;
            };

            void set_adof_owning_processor( const moris::sint aOwningProcessor )
            {
                mOwningProcessor = aOwningProcessor;
            };

            moris::moris_id get_adof_id()
            {
                return mAdofId;
            };

            moris::moris_id get_adof_external_id()
            {
                return mAdofExternalId;
            };

            moris::moris_index get_adof_external_ind()
            {
                return mAdofExternalInd;
            };

            moris::moris_id get_adof_owning_processor()
            {
                return mOwningProcessor;
            };

            moris::sint get_adof_type_time_identifier()
            {
                return mAdofTypeTimeIdentifier;
            };

            void set_adof_type_time_identifier( const moris::sint aAdofTypeTimeIdentifier)
            {
                mAdofTypeTimeIdentifier = aAdofTypeTimeIdentifier;
            };

        };
    }
}



#endif /* SRC_FEM_CL_ADOF_HPP_ */
