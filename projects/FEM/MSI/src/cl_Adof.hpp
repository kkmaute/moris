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

    public:
        Adof()
        {};

        ~Adof()
        {};

        void set_adof_id( const moris::uint aAdofId )
        {
            mAdofId = aAdofId;
        };

        const moris::uint get_adof_id()
        {
            return mAdofId;
        };

    };
    }
}



#endif /* SRC_FEM_CL_ADOF_HPP_ */
