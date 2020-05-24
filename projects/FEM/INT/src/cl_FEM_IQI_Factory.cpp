/*
 * cl_FEM_IQI_Factory.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include <memory>
#include "assert.hpp"
//FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_IQI_Volume.hpp"
#include "cl_FEM_IQI_Strain_Energy.hpp"
#include "cl_FEM_IQI_Stress.hpp"
#include "cl_FEM_IQI_Analytic.hpp"
#include "cl_FEM_IQI_Dof.hpp"
#include "cl_FEM_IQI_Property.hpp"
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"
#include "cl_FEM_IQI_H1_Error_Analytic.hpp"
#include "cl_FEM_IQI_H1_Semi_Error.hpp"
#include "cl_FEM_IQI_J_Integral.hpp"
#include "cl_FEM_IQI_K1_SENT.hpp"
#include "cl_FEM_IQI_Volume_Fraction.hpp"
#include "cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp"
namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        std::shared_ptr< IQI > IQI_Factory::create_IQI( IQI_Type aIQIType )
        {
            switch( aIQIType )
            {
                case IQI_Type::DOF :
                    return std::make_shared< IQI_Dof >();

                case IQI_Type::PROPERTY :
                    return std::make_shared< IQI_Property >();

                case IQI_Type::VOLUME :
                    return std::make_shared< IQI_Volume >();

                case IQI_Type::VOLUME_FRACTION :
                    return std::make_shared< IQI_Volume_Fraction >();

                case IQI_Type::STRAIN_ENERGY :
                    return std::make_shared< IQI_Strain_Energy >();

                case IQI_Type::STRESS :
                    return std::make_shared< IQI_Stress >();

                case IQI_Type::ANALYTIC :
                    return std::make_shared< IQI_Analytic >();

                case IQI_Type::L2_ERROR_ANALYTIC :
                    return std::make_shared< IQI_L2_Error_Analytic >();

                case IQI_Type::H1_ERROR_ANALYTIC :
                    return std::make_shared< IQI_H1_Error_Analytic >();

                case IQI_Type::H1_SEMI_ERROR :
                    return std::make_shared< IQI_H1_Semi_Error >();

                case IQI_Type::J_INTEGRAL :
                    return std::make_shared< IQI_J_Integral >();

                case IQI_Type::K1_SENT :
                    return std::make_shared< IQI_K1_SENT >();

                case IQI_Type::DRAG_COEFF :
                    return std::make_shared< IQI_Drag_Lift_Coefficient >( 1 );

                case IQI_Type::LIFT_COEFF :
                    return std::make_shared< IQI_Drag_Lift_Coefficient >( -1 );

                default:
                    MORIS_ERROR( false, " IQI_Factory::create_IQI - No IQI type specified. " );
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */



