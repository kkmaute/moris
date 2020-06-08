#ifndef MORIS_CL_GEN_INTEGRATION_PDV_HOST_HPP
#define MORIS_CL_GEN_INTEGRATION_PDV_HOST_HPP

#include "cl_GEN_Pdv_Enums.hpp"
#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {
        class GEN_Geometry_Object;

        class Integration_Pdv_Host
        {
        private :
            // Identifies the host node
            uint mNodeIndex;
            const Matrix<DDRMat>& mCoordinates;

            // PDV Information
            Cell<std::shared_ptr<Pdv>> mPdvs = Cell<std::shared_ptr<Pdv>>(3, nullptr);
            uint mStartingGlobalIndex;

        public:

            /**
             * Constructor
             *
             * @param aPdvTypes PDV types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV types
             */
            Integration_Pdv_Host(uint aNodeIndex,
                                 const Matrix<DDRMat>& aCoordinates,
                                 const Cell<PDV_Type>& aPdvTypes,
                                 uint aStartingGlobalIndex);

            /**
             * destructor
             */
            ~Integration_Pdv_Host();

            /**
             * Create PDV with real value
             *
             * @param aPdvType PDV type
             * @param aPdvVal PDV value
             */
            void create_pdv(PDV_Type aPdvType, moris::real aPdvVal);

            /**
             * Create PDV based on an intersection object
             *
             * @param aPdvType PDV type
             * @param aIntersection Pointer to an object with intersection information
             * @param aDimension 0, 1, or 2, corresponding to X, Y, and Z dimensions
             */
            void create_pdv(PDV_Type aPdvType, GEN_Geometry_Object* aIntersection, uint aDimension);

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
