#ifndef MORIS_CL_GEN_PDV_HOST_HPP_
#define MORIS_CL_GEN_PDV_HOST_HPP_

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
            Cell<std::shared_ptr<GEN_Pdv>> mPdvList;
            moris::map<PDV, uint> mPdvTypeMap;
            Matrix<DDUMat> mGlobalPdvIndices;
            Cell<bool> mActivePdvs;
            
            // list of properties for the dv types
            Cell< std::shared_ptr< GEN_Property > > mPdvProperties;
            
        public:
            
            /**
             * Constructor
             *
             * @param aPdvTypes PDV types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV types
             */
            Pdv_Host(const Cell<PDV>& aPdvTypes, uint aStartingGlobalIndex);
            
            /**
             * destructor
             */
            ~Pdv_Host();

            /**
             * Update the pdv type list to include potentially new PDV types
             *
             * @param aPdvTypes Potentially new PDV types to be added
             * @param aStartingGlobalIndex Global index to start assigning to PDV types
             * @return Number of added PDV types (unique PDVs)
             */
            uint add_pdv_types(const Cell<PDV>& aPdvTypes, uint aStartingGlobalIndex);
            
            /**
             * Create PDV with GEN field
             *
             * @param aPdvType PDV type
             * @param aFieldPointer Pointer to a GEN field
             * @param aNodeIndex Node index for pulling a value from the field
             */
            void create_pdv(PDV aPdvType, std::shared_ptr<GEN_Field> aFieldPointer, uint aNodeIndex);
            
            /**
             * Create PDV with GEN property
             *
             * @param aPdvType PDV type
             * @param aPropertyPointer Pointer to a GEN property
             */
            void create_pdv(PDV aPdvType, std::shared_ptr<GEN_Property> aPropertyPointer);
            
            /**
             * Create PDV with real value
             *
             * @param aPdvType PDV type
             * @param aPdvVal PDV value
             */
            void create_pdv(PDV aPdvType, moris::real aPdvVal);
            
            /**
             * Check if PDV type is active on this host
             *
             * @param aPdvType PDV type
             * @return if PDV type is active
             */
            bool is_active_type(PDV aPdvType);

            /**
             * Mark a pdv as being inactive
             *
             * @param aPdvType PDV type
             * @param aGlobalPdvTypeMap
             */
            void mark_pdv_as_inactive(PDV aPdvType);
            
            /**
             * Get global index for pdv by type
             *
             * @param aPdvType PDV type
             * @return Global index
             */
            uint get_global_index_for_pdv_type(PDV aPdvType);

            /**
             * Get the value of a PDV by type
             *
             * @param aPdvType PDV type
             * @return PDV value
             */
            real get_pdv_value(PDV aPdvType);
            
        };
    }  // end ge namepsace
} // end moris namespace

#endif /* MORIS_CL_GEN_PDV_HOST_HPP_ */
