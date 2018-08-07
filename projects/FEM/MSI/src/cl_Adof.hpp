/*
 * cl_Adof.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_ADOF_HPP_
#define SRC_FEM_CL_ADOF_HPP_

namespace moris
{
    namespace MSI
    {
    class Adof
    {
    private:
        moris::uint mAdofId;

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

        void set_adof_owning_processor( const moris::sint aOwningProcessor )
        {
            mOwningProcessor = aOwningProcessor;
        };

        const moris::uint get_adof_id()
        {
            return mAdofId;
        };

        const moris::uint get_adof_owning_processor()
        {
            return mOwningProcessor;
        };

    };
    }
}



#endif /* SRC_FEM_CL_ADOF_HPP_ */
