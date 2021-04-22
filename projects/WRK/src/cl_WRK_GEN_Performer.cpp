
#include "cl_WRK_GEN_Performer.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace wrk
    {

        Gen_Performer::Gen_Performer( std::shared_ptr<moris::ge::Geometry_Engine> aGeometryEngine )
                : mGeometryEngine(aGeometryEngine)
        {

        }

        uint 
        Gen_Performer::get_num_refinement_fields()
        {
            return mGeometryEngine->get_num_refinement_fields();
        } 

        real
        Gen_Performer::get_field_value(
                        uint                  aFieldIndex,
                        uint                  aNodeIndex,
                        const Matrix<DDRMat>& aCoordinates)
        {
            return mGeometryEngine->get_field_value(aFieldIndex,aNodeIndex,aCoordinates);
        }

        const Matrix< DDSMat > & 
        Gen_Performer::get_num_refinements( uint aFieldIndex )
        {
            return mGeometryEngine->get_num_refinements(aFieldIndex);
        }

        const Matrix< DDSMat > &
        Gen_Performer::get_refinement_mesh_indices(uint aFieldIndex )
        {
            return mGeometryEngine->get_refinement_mesh_indices(aFieldIndex);
        }

        sint
        Gen_Performer::get_refinement_function_index(
                        uint aFieldIndex,
                        uint aRefinementIndex)
        {
            return mGeometryEngine->get_refinement_function_index(aFieldIndex,aRefinementIndex);
        }
    }
}

