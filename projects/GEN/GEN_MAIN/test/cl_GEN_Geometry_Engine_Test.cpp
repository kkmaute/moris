#include "cl_GEN_Geometry_Engine_Test.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine_Test::Geometry_Engine_Test(
                Cell< std::shared_ptr<Geometry> > aGeometry,
                mtk::Interpolation_Mesh*          aMesh)
        : Geometry_Engine(aGeometry, aMesh)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Geometry> Geometry_Engine_Test::get_geometry(uint aGeometryIndex)
        {
            return mGeometries(aGeometryIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}