/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Library_IO_Standard.cpp
 *
 */

#include "cl_Library_IO_Standard.hpp"
#include "enums.hpp"
#include "parameters.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::Library_IO_Standard()
            : Library_IO()    // initialize base class data as usual
    {
    }

    //------------------------------------------------------------------------------------------------------------------

    Library_IO_Standard::~Library_IO_Standard()
    {
        // do nothing extra
    }

    //------------------------------------------------------------------------------------------------------------------

    bool Library_IO_Standard::is_module_supported( moris::Module_Type aModuleType )
    {
        switch ( aModuleType )
        {
            case Module_Type::OPT:
            case Module_Type::HMR:
            case Module_Type::STK:
            case Module_Type::XTK:
            case Module_Type::GEN:
            case Module_Type::FEM:
            case Module_Type::SOL:
            case Module_Type::MSI:
            case Module_Type::VIS:
            case Module_Type::MIG:
            case Module_Type::WRK:
            case Module_Type::MORISGENERAL:
                return true;
            default:
                return false;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
