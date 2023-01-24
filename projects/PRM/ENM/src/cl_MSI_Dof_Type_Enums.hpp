/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Dof_Type_Enums.hpp
 *
 */

#ifndef SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_
#define SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_

#include "cl_Map.hpp"

namespace moris
{
    namespace MSI
    {
        enum class Dof_Type
        {
            UX,    //< X-Displacement
            UY,    //< Y-Displacement
            UZ,    //< Z-Displacement

            TEMP,    //< Temperature

            L2,             //< Least Squares type
            MAPPING_DOF,    //< L2 Mapping
            LS1,            //< Level set
            LS2,            //< Level set

            NLSX,    //< X-Level set normal
            NLSY,    //< Y-Level set normal
            NLSZ,    //< Z-Level set normal

            THETA,    //< Heat method temperature
            PHID,     //< Heat method distance field
            PHISD,    //< Heat method signed distance field

            VX,    //< X-Velocity
            VY,    //< Y-Velocity
            VZ,    //< Z-Velocity

            P,      //< Pressure
            RHO,    //< Density

            VISCOSITY,     //< Turbulence viscosity
            STRESS_DOF,    //< Stress Dof

            E,     //< Energy density
            MX,    //< X-Momentum
            MY,    //< Y-Momentum
            MZ,    //< Z-Momentum

            EVP,    //< Entropy Variable for TD System Potential
            EVT,    //< Entropy production rate
            EVX,    //< X-Entropy transport of production rate
            EVY,    //< Y-Entropy transport of production rate
            EVZ,    //< Z-Entropy transport of production rate

            LAMBDAX,    //< x component of the Lagrange multiplier traction
            LAMBDAY,    //< y component of the Lagrange multiplier traction
            LAMBDAZ,    //< z component of the Lagrange multiplier traction

            UNDEFINED,    //< Undefined
            END_ENUM      //
        };

        //------------------------------------------------------------------------------

        inline map< std::string, enum MSI::Dof_Type >
        get_msi_dof_type_map()
        {
            map< std::string, enum MSI::Dof_Type > tMSIDofTypeMap;

            tMSIDofTypeMap[ "TEMP" ] = MSI::Dof_Type::TEMP;
            tMSIDofTypeMap[ "P" ]    = MSI::Dof_Type::P;
            tMSIDofTypeMap[ "RHO" ]  = MSI::Dof_Type::RHO;
            tMSIDofTypeMap[ "E" ]    = MSI::Dof_Type::E;
            tMSIDofTypeMap[ "EVP" ]  = MSI::Dof_Type::EVP;
            tMSIDofTypeMap[ "EVT" ]  = MSI::Dof_Type::EVT;

            tMSIDofTypeMap[ "UX" ]   = MSI::Dof_Type::UX;
            tMSIDofTypeMap[ "UY" ]   = MSI::Dof_Type::UY;
            tMSIDofTypeMap[ "UZ" ]   = MSI::Dof_Type::UZ;
            tMSIDofTypeMap[ "VX" ]   = MSI::Dof_Type::VX;
            tMSIDofTypeMap[ "VY" ]   = MSI::Dof_Type::VY;
            tMSIDofTypeMap[ "VZ" ]   = MSI::Dof_Type::VZ;
            tMSIDofTypeMap[ "MX" ]   = MSI::Dof_Type::MX;
            tMSIDofTypeMap[ "MY" ]   = MSI::Dof_Type::MY;
            tMSIDofTypeMap[ "MZ" ]   = MSI::Dof_Type::MZ;
            tMSIDofTypeMap[ "EVX" ]  = MSI::Dof_Type::EVX;
            tMSIDofTypeMap[ "EVY" ]  = MSI::Dof_Type::EVY;
            tMSIDofTypeMap[ "EVZ" ]  = MSI::Dof_Type::EVZ;
            tMSIDofTypeMap[ "NLSX" ] = MSI::Dof_Type::NLSX;
            tMSIDofTypeMap[ "NLSY" ] = MSI::Dof_Type::NLSY;
            tMSIDofTypeMap[ "NLSZ" ] = MSI::Dof_Type::NLSZ;

            tMSIDofTypeMap[ "L2" ]          = MSI::Dof_Type::L2;
            tMSIDofTypeMap[ "MAPPING_DOF" ] = MSI::Dof_Type::MAPPING_DOF;
            tMSIDofTypeMap[ "LS1" ]         = MSI::Dof_Type::LS1;
            tMSIDofTypeMap[ "LS2" ]         = MSI::Dof_Type::LS2;

            tMSIDofTypeMap[ "THETA" ] = MSI::Dof_Type::THETA;
            tMSIDofTypeMap[ "PHID" ]  = MSI::Dof_Type::PHID;
            tMSIDofTypeMap[ "PHISD" ] = MSI::Dof_Type::PHISD;

            tMSIDofTypeMap[ "VISCOSITY" ]  = MSI::Dof_Type::VISCOSITY;
            tMSIDofTypeMap[ "STRESS_DOF" ] = MSI::Dof_Type::STRESS_DOF;

            tMSIDofTypeMap[ "UNDEFINED" ] = MSI::Dof_Type::UNDEFINED;

            return tMSIDofTypeMap;
        }

        //------------------------------------------------------------------------------

        inline map< enum MSI::Dof_Type, std::string >
        get_dof_type_name_map()
        {
            map< enum MSI::Dof_Type, std::string > tMSIDofTypeMap;

            tMSIDofTypeMap[ MSI::Dof_Type::TEMP ] = "TEMP";
            tMSIDofTypeMap[ MSI::Dof_Type::P ]    = "P";
            tMSIDofTypeMap[ MSI::Dof_Type::RHO ]  = "RHO";
            tMSIDofTypeMap[ MSI::Dof_Type::RHO ]  = "RHO";
            tMSIDofTypeMap[ MSI::Dof_Type::E ]    = "E";
            tMSIDofTypeMap[ MSI::Dof_Type::EVP ]  = "EVP";
            tMSIDofTypeMap[ MSI::Dof_Type::EVT ]  = "EVT";

            tMSIDofTypeMap[ MSI::Dof_Type::UX ]   = "UX";
            tMSIDofTypeMap[ MSI::Dof_Type::UY ]   = "UY";
            tMSIDofTypeMap[ MSI::Dof_Type::UZ ]   = "UZ";
            tMSIDofTypeMap[ MSI::Dof_Type::VX ]   = "VX";
            tMSIDofTypeMap[ MSI::Dof_Type::VY ]   = "VY";
            tMSIDofTypeMap[ MSI::Dof_Type::VZ ]   = "VZ";
            tMSIDofTypeMap[ MSI::Dof_Type::MX ]   = "MX";
            tMSIDofTypeMap[ MSI::Dof_Type::MY ]   = "MY";
            tMSIDofTypeMap[ MSI::Dof_Type::MZ ]   = "MZ";
            tMSIDofTypeMap[ MSI::Dof_Type::EVX ]  = "EVX";
            tMSIDofTypeMap[ MSI::Dof_Type::EVY ]  = "EVY";
            tMSIDofTypeMap[ MSI::Dof_Type::EVZ ]  = "EVZ";
            tMSIDofTypeMap[ MSI::Dof_Type::NLSX ] = "NLSX";
            tMSIDofTypeMap[ MSI::Dof_Type::NLSY ] = "NLSY";
            tMSIDofTypeMap[ MSI::Dof_Type::NLSZ ] = "NLSZ";

            tMSIDofTypeMap[ MSI::Dof_Type::L2 ]          = "L2";
            tMSIDofTypeMap[ MSI::Dof_Type::MAPPING_DOF ] = "MAPPING_DOF";
            tMSIDofTypeMap[ MSI::Dof_Type::LS1 ]         = "LS1";
            tMSIDofTypeMap[ MSI::Dof_Type::LS2 ]         = "LS2";

            tMSIDofTypeMap[ MSI::Dof_Type::THETA ] = "THETA";
            tMSIDofTypeMap[ MSI::Dof_Type::PHID ]  = "PHID";
            tMSIDofTypeMap[ MSI::Dof_Type::PHISD ] = "PHISD";

            tMSIDofTypeMap[ MSI::Dof_Type::VISCOSITY ]  = "VISCOSITY";
            tMSIDofTypeMap[ MSI::Dof_Type::STRESS_DOF ] = "STRESS_DOF";

            tMSIDofTypeMap[ MSI::Dof_Type::UNDEFINED ] = "UNDEFINED";

            return tMSIDofTypeMap;
        }
    }    // namespace MSI
}    // namespace moris

#endif /* SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_ */
