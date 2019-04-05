/*
 * cl_GE_Factory.cpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */
#include "cl_GE_Factory.hpp"


namespace moris
{
    namespace ge
    {
        Ge_Factory::Ge_Factory(){}
        Ge_Factory::~Ge_Factory(){}

        std::shared_ptr< Geometry > Ge_Factory::set_geometry_type( const enum type aGeomType )
        {
            std::shared_ptr< Geometry > tGeomPointer = nullptr;
            switch(aGeomType)
            {
            case(type::ANALYTIC):
                    tGeomPointer = std::make_shared< Analytical >();
                    break;
            case(type::DISCRETE):
                    tGeomPointer = std::make_shared< Discrete >();
                    break;
            default:
                    MORIS_ERROR(false, "Ge_Factory::set_geometry_type() please input a valid geometry type");
                    break;
            }
        return tGeomPointer;
        }
    } /* namespace gen */
} /* namespace moris */
