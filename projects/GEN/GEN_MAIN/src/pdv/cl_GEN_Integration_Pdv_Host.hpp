#ifndef MORIS_CL_GEN_INTEGRATION_PDV_HOST_HPP
#define MORIS_CL_GEN_INTEGRATION_PDV_HOST_HPP

#include "cl_GEN_Pdv_Enums.hpp"
#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {
        // Forward declaration for struct
        class Geometry;
        class Integration_Pdv_Host;

        // Intersection data structure
        struct Intersection
        {
            std::shared_ptr<Geometry> mGeometry;
            Cell<std::shared_ptr<Integration_Pdv_Host>> mNodes;
        };

        class Integration_Pdv_Host
        {
        private :
            // Identifies the host node
            uint mNodeIndex;
            const Matrix<DDRMat>& mCoordinates;

            // PDV Information
            uint mStartingGlobalIndex;
            std::shared_ptr<Intersection> mIntersection = nullptr;

        public:

            /**
             * Constructor
             *
             * @param aPdvTypes PDV types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV types
             */
            Integration_Pdv_Host(uint aNodeIndex,
                                 const Matrix<DDRMat>& aCoordinates,
                                 uint aStartingGlobalIndex,
                                 std::shared_ptr<Intersection> aIntersection = nullptr);

            /**
             * destructor
             */
            ~Integration_Pdv_Host();

            /**
             * Sets this PDV host as an intersection point
             *
             * @param aIntersection Information about other PDV hosts this PDV hosts depends on based on an intersection
             */
            void set_as_intersection(std::shared_ptr<Intersection> aIntersection);

            /**
             * Check if PDV type is active on this host
             *
             * @param aPdvType PDV type
             * @return if PDV type is active
             */
            bool is_active_type(PDV_Type aPdvType);

            /**
             * Get global index for pdv by type
             *
             * @param aPdvType PDV type
             * @return Global index
             */
            uint get_global_index_for_pdv_type(PDV_Type aPdvType);

            /**
             * Get the starting global index for the integration PDVs
             *
             * @return The global index of the first PDV on the host
             */
            uint get_starting_global_index();

            /**
             * Get the value of a PDV by type
             *
             * @param aPdvType PDV type
             * @return Value on this PDV
             */
            real get_pdv_value(PDV_Type aPdvType);

            /**
             * Gets all of the sensitivity vectors on each PDV
             *
             * @param aSensitivities
             */
            void get_all_sensitivities(Matrix<DDRMat>& aSensitivities);

        };
    }
}

#endif //MORIS_CL_GEN_INTEGRATION_PDV_HOST_HPP
