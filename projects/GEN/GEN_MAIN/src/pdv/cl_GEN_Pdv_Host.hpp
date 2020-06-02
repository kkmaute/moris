#ifndef MORIS_CL_GEN_PDV_Type_HOST_HPP_
#define MORIS_CL_GEN_PDV_Type_HOST_HPP_

// GEN_MAIN
#include "cl_GEN_Pdv.hpp"
#include "cl_GEN_Pdv_Property.hpp"

// GEN_CORE
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Host
        {
        private :
            // Identifies the host node
            uint mNodeIndex;
            Matrix<DDRMat> mCoordinates;

            // Information about the contained PDVs
            Cell<std::shared_ptr<Pdv>> mPdvList;
            moris::map<PDV_Type, uint> mPdvTypeMap;
            Matrix<DDUMat> mGlobalPdvIndices;
            Cell<bool> mActivePdvs;
            
        public:
            
            /**
             * Constructor
             *
             * @param aPdvTypes PDV_Type types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV_Type types
             */
            Pdv_Host(uint aNodeIndex, const Cell<PDV_Type>& aPdvTypes, uint aStartingGlobalIndex);
            
            /**
             * destructor
             */
            ~Pdv_Host();

            /**
             * Update the pdv type list to include potentially new PDV_Type types
             *
             * @param aPdvTypes Potentially new PDV_Type types to be added
             * @param aStartingGlobalIndex Global index to start assigning to PDV_Type types
             * @return Number of added PDV_Type types (unique PDVs)
             */
            uint add_pdv_types(const Cell<PDV_Type>& aPdvTypes, uint aStartingGlobalIndex);
            
            /**
             * Create PDV_Type with GEN property
             *
             * @param aPdvType PDV_Type type
             * @param aPropertyPointer Pointer to a GEN property
             */
            void create_pdv(PDV_Type aPdvType, std::shared_ptr<Property> aPropertyPointer);
            
            /**
             * Create PDV_Type with real value
             *
             * @param aPdvType PDV_Type type
             * @param aPdvVal PDV_Type value
             */
            void create_pdv(PDV_Type aPdvType, moris::real aPdvVal);
            
            /**
             * Check if PDV_Type type is active on this host
             *
             * @param aPdvType PDV_Type type
             * @return if PDV_Type type is active
             */
            bool is_active_type(PDV_Type aPdvType);

            /**
             * Mark a pdv as being inactive
             *
             * @param aPdvType PDV_Type type
             * @param aGlobalPdvTypeMap
             */
            void mark_pdv_as_inactive(PDV_Type aPdvType);
            
            /**
             * Get global index for pdv by type
             *
             * @param aPdvType PDV_Type type
             * @return Global index
             */
            uint get_global_index_for_pdv_type(PDV_Type aPdvType);

            /**
             * Get all of the global PDV indices on this host
             *
             * @param aGlobalPdvIndices matrix of indices to be returned
             */
            const Matrix<DDUMat>& get_all_global_indices();

            /**
             * Get the value of a PDV_Type by type
             *
             * @param aPdvType PDV_Type type
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
    }  // end ge namepsace
} // end moris namespace

#endif /* MORIS_CL_GEN_PDV_HOST_HPP_ */
