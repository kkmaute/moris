#ifndef MORIS_CL_GEN_PDV_Type_HOST_HPP_
#define MORIS_CL_GEN_PDV_Type_HOST_HPP_

#include "cl_GEN_Pdv_Enums.hpp"
#include "cl_GEN_Pdv.hpp"

namespace moris
{
    namespace ge
    {
        class GEN_Geometry_Object;
        class Property;

        class Pdv_Host
        {
        private :
            // Identifies the host node
            uint mNodeIndex;
            Matrix<DDRMat> mCoordinates;

            // Information about the contained PDVs
            Cell<std::shared_ptr<Pdv>> mPdvs;
            moris::map<PDV_Type, uint> mPdvTypeMap;
            Matrix<DDUMat> mGlobalPdvIndices;
            
        public:
            
            /**
             * Constructor
             *
             * @param aPdvTypes PDV types for this host
             * @param aStartingGlobalIndex Global index to start assigning new PDV types
             */
            Pdv_Host(uint aNodeIndex, const Cell<PDV_Type>& aPdvTypes, uint aStartingGlobalIndex);
            
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
            uint add_pdv_types(const Cell<PDV_Type>& aPdvTypes, uint aStartingGlobalIndex);
            
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
             * Get all of the global PDV indices on this host
             *
             * @param aGlobalPdvIndices matrix of indices to be returned
             */
            const Matrix<DDUMat>& get_all_global_indices();

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
    }  // end ge namepsace
} // end moris namespace

#endif /* MORIS_CL_GEN_PDV_HOST_HPP_ */
