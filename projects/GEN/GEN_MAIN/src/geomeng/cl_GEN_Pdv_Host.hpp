#ifndef MORIS_CL_GEN_PDV_Type_HOST_HPP_
#define MORIS_CL_GEN_PDV_Type_HOST_HPP_

// GEN_MAIN
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Pdv.hpp"
#include "cl_GEN_Property.hpp"

// GEN_CORE
#include "cl_GEN_Pdv_Enums.hpp"

namespace moris
{
    namespace ge
    {
        class Pdv_Host
        {
            
        private :
            // Information about the contained PDVs
            Cell<std::shared_ptr<Pdv>> mPdvList;
            moris::map<PDV_Type, uint> mPdvTypeMap;
            Matrix<DDUMat> mGlobalPdvIndices;
            Cell<bool> mActivePdvs;
            
            // Coordinates
            Matrix<DDRMat> mCoordinates = {{0.0, 0.0, 0.0}};
            
        public:
            
            /**
             * Constructor
             *
             * @param aPdvTypes PDV_Type types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV_Type types
             */
            Pdv_Host(const Cell<PDV_Type>& aPdvTypes, uint aStartingGlobalIndex);
            
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
             * Create PDV_Type with GEN field
             *
             * @param aPdvType PDV_Type type
             * @param aFieldPointer Pointer to a GEN field
             * @param aNodeIndex Node index for pulling a value from the field
             */
            void create_pdv(PDV_Type aPdvType, std::shared_ptr<GEN_Field> aFieldPointer, uint aNodeIndex);
            
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
             * Get the value of a PDV_Type by type
             *
             * @param aPdvType PDV_Type type
             * @return PDV_Type value
             */
            real get_pdv_value(PDV_Type aPdvType);
            
        };
    }  // end ge namepsace
} // end moris namespace

#endif /* MORIS_CL_GEN_PDV_HOST_HPP_ */
