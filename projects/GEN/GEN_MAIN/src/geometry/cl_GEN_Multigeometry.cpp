#include "cl_GEN_Multigeometry.hpp"

namespace moris
{
    namespace ge
    {

        //--------------------------------------------------------------------------------------------------------------

        Multigeometry::Multigeometry(Cell<std::shared_ptr<Geometry>> aGeometries,
                                     std::string aID)
                : Field(Matrix<DDRMat>(0, 0),
                        aGeometries(0)->get_num_refinements(),
                        aGeometries(0)->get_refinement_function_index(),
                        aGeometries(0)->get_bspline_mesh_index(),
                        aGeometries(0)->get_bspline_lower_bound(),
                        aGeometries(0)->get_bspline_upper_bound()),
                  mGeometries(aGeometries),
                  mID(aID)
        {
            MORIS_ERROR(mGeometries.size() > 0, "A GEN Multigeometry must be created with at least one geometry.");
        }

        //--------------------------------------------------------------------------------------------------------------

        real Multigeometry::evaluate_field_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates)
        {
            real tResult = mGeometries(0)->evaluate_field_value(aNodeIndex, aCoordinates);
            for (uint tGeometryIndex = 1; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                tResult = std::min(tResult, mGeometries(tGeometryIndex)->evaluate_field_value(aNodeIndex, aCoordinates));
            }
            return tResult;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Multigeometry::evaluate_sensitivity(uint                  aNodeIndex,
                                                 const Matrix<DDRMat>& aCoordinates,
                                                 Matrix<DDRMat>&       aSensitivities)
        {
            // Find which geometry is the minimum
            real tMin = mGeometries(0)->evaluate_field_value(aNodeIndex, aCoordinates);
            uint tMinGeometryIndex = 0;
            for (uint tGeometryIndex = 1; tGeometryIndex < mGeometries.size(); tGeometryIndex++)
            {
                real tResult = mGeometries(tGeometryIndex)->evaluate_field_value(aNodeIndex, aCoordinates);
                if (tResult < tMin)
                {
                    tMin = tResult;
                    tMinGeometryIndex = tGeometryIndex;
                }
            }

            // Return relevant sensitivity
            mGeometries(tMinGeometryIndex)->evaluate_sensitivity(aNodeIndex, aCoordinates, aSensitivities);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Multigeometry::add_geometry(std::shared_ptr<Geometry> aGeometry)
        {
            mGeometries.push_back(aGeometry);
        }

        //--------------------------------------------------------------------------------------------------------------

        std::string Multigeometry::get_id()
        {
            return mID;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Multigeometry::evaluate_all_sensitivities(uint                  aNodeIndex,
                                                       const Matrix<DDRMat>& aCoordinates,
                                                       Matrix<DDRMat>&       aSensitivities)
        {
        }

        //--------------------------------------------------------------------------------------------------------------

    }
}
