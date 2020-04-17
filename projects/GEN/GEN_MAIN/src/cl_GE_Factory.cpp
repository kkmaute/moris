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

        std::shared_ptr< Geometry_Analytic > Ge_Factory::set_geometry_type(const enum GeomType aGeomType )
        {
            std::shared_ptr< Geometry_Analytic > tGeomPointer = nullptr;
            switch(aGeomType)
            {
            case(GeomType::ANALYTIC):
                    tGeomPointer = std::make_shared< Analytic >();
                    break;
            case(GeomType::DISCRETE):
                    tGeomPointer = std::make_shared< Discrete >();
                    break;
            case(GeomType::SDF):
                    tGeomPointer = std::make_shared< SDF >();
                    break;
            default:
                    MORIS_ERROR(false, "Ge_Factory::set_geometry_type() please input a valid geometry type");
                    break;
            }
        return tGeomPointer;
        }
    } /* namespace gen */
} /* namespace moris */
