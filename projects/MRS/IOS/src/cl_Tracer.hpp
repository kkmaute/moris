#ifndef MORIS_IOS_CL_TRACER_HPP_
#define MORIS_IOS_CL_TRACER_HPP_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include "cl_Logger.hpp"

//get access to the global clock in gLogger;
extern moris::Logger gLogger;

namespace moris
{
class Tracer
    {

    //-------------------------------- PRIVATE --------------------------------//
    private:

    // class has no member variables / is empty


    //-------------------------------- PUBLIC ---------------------------------//
    public:

        /**
         * Constructor
         *
         * @param aEntityBase Entity base
         * @param aEntityType Entity type
         * @param aEntityAction Entity action
         */
        Tracer(std::string aEntityBase,
               std::string aEntityType,
               std::string aEntityAction)
        {
            gLogger.sign_in( aEntityBase, aEntityType, aEntityAction );
        }

        /**
         * Constructor
         *
         * @param aEntityBase Entity base
         * @param aEntityAction Entity action
         */
        Tracer(std::string aEntityBase,
               std::string aEntityAction)
        {
            gLogger.sign_in( aEntityBase, "NoType", aEntityAction );
        }

        // destructor: automatically perform sign out operation when tracer gets destructed
        ~Tracer()
        {
            gLogger.sign_out();
        };


    }; // class Tracer
} // namespace moris

#endif  /* MORIS_IOS_CL_TRACER_HPP_ */

