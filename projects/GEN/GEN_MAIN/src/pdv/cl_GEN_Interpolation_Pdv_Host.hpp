#ifndef MORIS_CL_GEN_INTERPOLATION_PDV_Type_HOST_HPP_
#define MORIS_CL_GEN_INTERPOLATION_PDV_Type_HOST_HPP_

#include "cl_GEN_Pdv_Enums.hpp"
#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {
        class Property;
        class Pdv_Host_Manager;

        class Interpolation_Pdv_Host
        {
        private :
            Pdv_Host_Manager * mPdvHostManager;

            // Identifies the host node
            moris_index mNodeIndex;
            moris_id mNodeId;
            moris_index mNodeOwner;

            Matrix<DDRMat> mCoordinates;

            // Information about the contained PDVs
            Cell<std::shared_ptr<Pdv>> mPdvs;
            
        public:
            
            /**
             * Constructor
             *
             * @param aPdvTypes PDV types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV types
             */
            Interpolation_Pdv_Host(
                    Pdv_Host_Manager * aPdvHostManager,
                    const moris_index & aNodeIndex,
                    const moris_id & aNodeId,
                    const moris_index & aNodeOwner,
                    const Matrix<DDRMat> & aCoordinates,
                    const Cell<PDV_Type> & aPdvTypes);
            
            /**
             * destructor
             */
            ~Interpolation_Pdv_Host();
            
            moris_id get_pdv_vertex_id()
            {
                return mNodeId;
            };

            moris_index get_pdv_owning_processor()
            {
                return mNodeOwner;
            };

            void set_pdv_id(
                    enum PDV_Type aPdvType,
                    const moris_id aCounterId );

            moris_id get_pdv_id(
                    enum PDV_Type aPdvType);

            /**
             * Create PDV with real value
             *
             * @param aPdvType PDV type
             * @param aPdvVal PDV value
             */
            void create_pdv(PDV_Type aPdvType, moris::real aPdvVal);

            /**
             * Create PDV with GEN property
             *
             * @param aPdvType PDV type
             * @param aPropertyPointer Pointer to a GEN property
             */
            void create_pdv(PDV_Type aPdvType, std::shared_ptr<Property> aPropertyPointer);
            
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
            void set_global_index_for_pdv_type(PDV_Type aPdvType,  moris_id aId);

            /**
             * Get global index for pdv by type
             *
             * @param aPdvType PDV type
             * @return Global index
             */
            uint get_global_index_for_pdv_type(PDV_Type aPdvType);

            /**
             * Get all of the global PDV indices on this host
             *
             * @param aGlobalPdvIndices matrix of indices to be returned
             */
            Matrix<DDUMat> get_all_global_indices();

            /**
             * Get the value of a PDV by type
             *
             * @param aPdvType PDV type
             * @return Value on this PDV
             */
            real get_pdv_value(PDV_Type aPdvType);

            bool get_pdv_exists(PDV_Type aPdvType);

            /**
             * Gets all of the sensitivity vectors on each PDV
             *
             * @param aSensitivities
             */
            void get_all_sensitivities(Matrix<DDRMat>& aSensitivities);
            
        };
    }
}

#endif /* MORIS_CL_GEN_INTERPOLATION_PDV_HOST_HPP_ */
