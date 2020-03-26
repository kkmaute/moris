/*
 * cl_FEM_IQI_Factory.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include <memory>
#include "assert.hpp"
#include "cl_FEM_IQI_Factory.hpp"       //FEM/INT/src
#include "cl_FEM_IQI_Volume.hpp"        //FEM/INT/src
#include "cl_FEM_IQI_Strain_Energy.hpp" //FEM/INT/src
#include "cl_FEM_IQI_Stress.hpp"        //FEM/INT/src
#include "cl_FEM_IQI_Analytic.hpp"
#include "cl_FEM_IQI_Dof.hpp"           //FEM/INT/src
#include "cl_FEM_IQI_Property.hpp"           //FEM/INT/src
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"           //FEM/INT/src
#include "cl_FEM_IQI_H1_Error_Analytic.hpp"           //FEM/INT/src
#include "cl_FEM_IQI_H1_Semi_Error.hpp"           //FEM/INT/src
#include "cl_FEM_IQI_J_Integral.hpp"           //FEM/INT/src
#include "cl_FEM_IQI_Volume_Fraction.hpp"           //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        std::shared_ptr< IQI > IQI_Factory::create_IQI( IQI_Type aIQIType )
        {
            switch( aIQIType )
            {
                case ( IQI_Type::DOF ):
                    return std::make_shared< IQI_Dof >();

                case ( IQI_Type::PROPERTY ):
                    return std::make_shared< IQI_Property >();

                case ( IQI_Type::VOLUME ):
                    return std::make_shared< IQI_Volume >();

                case ( IQI_Type::VOLUME_FRACTION ):
                    return std::make_shared< IQI_Volume_Fraction >();

                case ( IQI_Type::STRAIN_ENERGY ):
                    return std::make_shared< IQI_Strain_Energy >();
                case( IQI_Type::STRESS):
                    return std::make_shared< IQI_Stress >();
                case( IQI_Type::ANALYTIC):
                    return std::make_shared< IQI_Analytic >();

                case ( IQI_Type::L2_ERROR_ANALYTIC ):
                    return std::make_shared< IQI_L2_Error_Analytic >();

                case ( IQI_Type::H1_ERROR_ANALYTIC ):
                    return std::make_shared< IQI_H1_Error_Analytic >();

                case ( IQI_Type::H1_SEMI_ERROR ):
                    return std::make_shared< IQI_H1_Semi_Error >();

                case ( IQI_Type::J_INTEGRAL ):
                    return std::make_shared< IQI_J_Integral >();

                default:
                    MORIS_ERROR( false, " IQI_Factory::create_IQI - No IQI type specified. " );
                    return nullptr;
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */



