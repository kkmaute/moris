#ifndef MORIS_CL_GEN_PDV_INTERSECTION_HPP
#define MORIS_CL_GEN_PDV_INTERSECTION_HPP

#include "cl_GEN_Pdv.hpp"
#include "cl_GEN_Geometry_Object.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Intersection : public Pdv
        {

        private:
            GEN_Geometry_Object* mIntersection;
            uint mDimension;

        public:
            /**
             * Constructor
             *
             * @param aPropertyPointer a GEN property pointer
             */
            Pdv_Intersection(GEN_Geometry_Object* aIntersection, uint aDimension);

            /**
             * Get the PDV value
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @return Current value of this PDV
             */
            real get_value(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates);

            /**
             * Get the PDV sensitivity with respect to ADVs
             *
             * @param aNodeIndex Node index
             * @param aCoordinates Coordinate values
             * @param aSensitivities Matrix of sensitivities to be returned
             */
            void get_sensitivity(uint aNodeIndex, const Matrix<DDRMat>& aCoordinates, Matrix<DDRMat>& aSensitivities);

        };
    }
}


#endif //MORIS_CL_GEN_PDV_INTERSECTION_HPP
