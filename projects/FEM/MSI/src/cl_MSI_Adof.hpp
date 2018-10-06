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
        moris::uint mAdofId = -1;
        moris::uint mAdofExternalId; //FIXME delte
        moris::uint mAdofExternalInd; //external Ind. only for HMR use

        moris::uint mOwningProcessor;

    public:
        Adof()
        {};

        ~Adof()
        {};

        void set_adof_id( const moris::uint aAdofId )
        {
            mAdofId = aAdofId;
        };

        void set_adof_external_id( const moris::uint aAdofExtId )   //FIXME delte
        {
            mAdofExternalId = aAdofExtId;
        };

        void set_adof_external_ind( const moris::uint aAdofExtInd )
        {
            mAdofExternalInd = aAdofExtInd;
        };

        void set_adof_owning_processor( const moris::sint aOwningProcessor )
        {
            mOwningProcessor = aOwningProcessor;
        };

        moris::uint get_adof_id()
        {
            return mAdofId;
        };

        moris::uint get_adof_external_id()
        {
            return mAdofExternalId;
        };

        moris::uint get_adof_external_ind()
        {
            return mAdofExternalInd;
        };

        moris::moris_id get_adof_owning_processor()
        {
            return mOwningProcessor;
        };

    };
    }
}



#endif /* SRC_FEM_CL_ADOF_HPP_ */
