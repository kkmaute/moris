#include "cl_GEN_Pdv_Host.hpp"

namespace moris
{
    namespace ge
    {
        
        //--------------------------------------------------------------------------------------------------------------
        
        Pdv_Host::Pdv_Host(const Cell<PDV_Type>& aPdvTypes, uint aGlobalIndex)
        : mPdvList(aPdvTypes.size(), nullptr),
          mGlobalPdvIndices(aPdvTypes.size(), 1),
          mActivePdvs(aPdvTypes.size(), true)
        {
            for (uint tPdvIndex = 0; tPdvIndex < aPdvTypes.size(); tPdvIndex++)
            {
                // Map from PDV_Type type to local index
                mPdvTypeMap[aPdvTypes(tPdvIndex)] = tPdvIndex;

                // local index to global index
                mGlobalPdvIndices(tPdvIndex) = aGlobalIndex++;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
        Pdv_Host::~Pdv_Host()
        {
        }

        //--------------------------------------------------------------------------------------------------------------
        
        uint Pdv_Host::add_pdv_types(const Cell<PDV_Type>& aPdvTypes, uint aGlobalIndex)
        {
            // Check for existing PDVs and add new ones to map
            uint tOriginalPdvs = mPdvTypeMap.size();
            uint tPdvIndex = tOriginalPdvs;
            for (uint tNewPdvIndex = 0; tNewPdvIndex < aPdvTypes.size(); tNewPdvIndex++)
            {
                if (!mPdvTypeMap.key_exists(aPdvTypes(tNewPdvIndex)))
                {
                    mPdvTypeMap[aPdvTypes(tNewPdvIndex)] = tPdvIndex;
                    tPdvIndex++;
                }
            }

            // Extend PDV_Type list and indices
            mPdvList.resize(tPdvIndex);
            mGlobalPdvIndices.resize(tPdvIndex, 1);
            mActivePdvs.resize(tPdvIndex, true);
            for (uint tNewPdvIndex = 0; tNewPdvIndex < tPdvIndex - tOriginalPdvs; tNewPdvIndex++)
            {
                mGlobalPdvIndices(tOriginalPdvs + tNewPdvIndex) = aGlobalIndex++;
            }

            // Return number of added PDV_Type types
            return tPdvIndex - tOriginalPdvs;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Pdv_Host::create_pdv(PDV_Type aPdvType, std::shared_ptr<GEN_Field> aFieldPointer, uint aNodeIndex)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.create_pdv() using GEN field with PDV_Type type that doesn't exist on this host.");

            // create a pdv with field pointer
            mPdvList(mPdvTypeMap[aPdvType]) = std::make_shared<GEN_Pdv>(aFieldPointer, aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host::create_pdv(PDV_Type aPdvType, std::shared_ptr<GEN_Property> aPropertyPointer)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.create_pdv() using GEN property with PDV_Type type that doesn't exist on this host.");

            // create a pdv with property pointer
            mPdvList(mPdvTypeMap[aPdvType]) = std::make_shared< GEN_Pdv >( aPropertyPointer );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host::create_pdv(PDV_Type aPdvType, moris::real aPdvVal)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.create_pdv() using pdv value with PDV_Type type that doesn't exist on this host.");

            // create a pdv with pdv value
            mPdvList(mPdvTypeMap[aPdvType]) = std::make_shared< GEN_Pdv >( aPdvVal );
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Pdv_Host::is_active_type(PDV_Type aPdvType)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.is_active_type() with type that doesn't exist on this host.");

            // return if active PDV_Type
            return mActivePdvs(mPdvTypeMap[aPdvType]);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host::mark_pdv_as_inactive(PDV_Type aPdvType)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.mark_pdv_as_inactive() with type that doesn't exist on this host.");

            // Set PDV_Type to be inactive
            mActivePdvs(mPdvTypeMap[aPdvType]) = false;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Pdv_Host::get_global_index_for_pdv_type(PDV_Type aPdvType)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                         "Tried to call Pdv_Host.get_global_index_for_pdv_type() with type that doesn't exist on this host.");

            // Return id from map
            return mGlobalPdvIndices(mPdvTypeMap[aPdvType]);
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv_Host::get_pdv_value(PDV_Type aPdvType)
        {
            // Check PDV_Type type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                         "Tried to call Pdv_Host.get_pdv_value() with type that doesn't exist on this host.");

            // Return value
            return mPdvList(mPdvTypeMap[aPdvType])->get_val()(0);
        }

        //--------------------------------------------------------------------------------------------------------------
        
    } // end ge namepsace
} // end moris namespace

