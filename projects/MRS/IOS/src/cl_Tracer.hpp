/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Tracer.hpp
 *
 */

#ifndef MORIS_IOS_CL_TRACER_HPP_
#define MORIS_IOS_CL_TRACER_HPP_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cstring>

#include "cl_Logger.hpp"
#include "Log_Constants.hpp"

//get access to the global clock in gLogger;
extern moris::Logger gLogger;

namespace moris
{
class Tracer
    {

    //-------------------------------- PRIVATE --------------------------------//
    private:

    // class has no member variables / is empty
    bool mSignOut = true;

    //-------------------------------- PUBLIC ---------------------------------//
    public:

        /**
         * Constructor for tracer
         *
         * @param aEntityBase Entity base
         * @param aEntityType Entity type
         * @param aEntityAction Entity action
         */
        Tracer(std::string aEntityBase,
               std::string aEntityType,
               std::string aEntityAction,
               moris::uint aThreshold = 0,
               moris::uint aLevel = 0)
        {
            if(aLevel <= aThreshold)
            {
                gLogger.sign_in( aEntityBase, aEntityType, aEntityAction );
            }
            else
            {
                mSignOut = false;
            }
        }

        /**
         * Constructor
         * Use if no entity type can be specified
         *
         * @param aEntityBase Entity base
         * @param aEntityAction Entity action
         */
        Tracer(std::string aEntityBase,
               std::string aEntityAction)
        {
            gLogger.sign_in( aEntityBase, LOGGER_NON_SPECIFIC_ENTITY_TYPE, aEntityAction );
        }

        // destructor: automatically perform sign out operation when tracer gets destructed
        ~Tracer()
        {
            if(mSignOut){gLogger.sign_out();}
        };

    }; // class Tracer
} // namespace moris

#endif  /* MORIS_IOS_CL_TRACER_HPP_ */

