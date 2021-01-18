#include "cl_GEN_Geometry_Engine_Test.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Geometry_Engine_Test::Geometry_Engine_Test(
                mtk::Interpolation_Mesh*   aMesh,
                Geometry_Engine_Parameters aParameters)
                : Geometry_Engine(aMesh, aParameters)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Geometry> Geometry_Engine_Test::get_geometry(uint aGeometryIndex)
        {
            return mGeometries(aGeometryIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        std::shared_ptr<Property> Geometry_Engine_Test::get_property(uint aPropertyIndex)
        {
            return mProperties(aPropertyIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}