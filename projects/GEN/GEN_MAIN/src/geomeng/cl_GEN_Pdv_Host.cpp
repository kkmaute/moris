#include "cl_GEN_Pdv_Host.hpp"

namespace moris
{
    namespace ge
    {
        
        //--------------------------------------------------------------------------------------------------------------
        
        Pdv_Host::Pdv_Host(const Cell<PDV>& aPdvTypes, uint aGlobalIndex)
        : mPdvList(aPdvTypes.size(), nullptr),
          mGlobalPdvIndices(aPdvTypes.size(), 1),
          mActivePdvs(aPdvTypes.size(), true)
        {
            for (uint tPdvIndex = 0; tPdvIndex < aPdvTypes.size(); tPdvIndex++)
            {
                // Map from PDV type to local index
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
        
        uint Pdv_Host::add_pdv_types(const Cell<PDV>& aPdvTypes, uint aGlobalIndex)
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

            // Extend PDV list and indices
            mPdvList.resize(tPdvIndex);
            mGlobalPdvIndices.resize(tPdvIndex, 1);
            mActivePdvs.resize(tPdvIndex, true);
            for (uint tNewPdvIndex = 0; tNewPdvIndex < tPdvIndex - tOriginalPdvs; tNewPdvIndex++)
            {
                mGlobalPdvIndices(tOriginalPdvs + tNewPdvIndex) = aGlobalIndex++;
            }

            // Return number of added PDV types
            return tPdvIndex - tOriginalPdvs;
        }

        //--------------------------------------------------------------------------------------------------------------
        
        void Pdv_Host::create_pdv(PDV aPdvType, std::shared_ptr<GEN_Field> aFieldPointer, uint aNodeIndex)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.create_pdv() using GEN field with PDV type that doesn't exist on this host.");

            // create a pdv with field pointer
            mPdvList(mPdvTypeMap[aPdvType]) = std::make_shared<GEN_Pdv>(aFieldPointer, aNodeIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host::create_pdv(PDV aPdvType, std::shared_ptr<GEN_Property> aPropertyPointer)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.create_pdv() using GEN property with PDV type that doesn't exist on this host.");

            // create a pdv with property pointer
            mPdvList(mPdvTypeMap[aPdvType]) = std::make_shared< GEN_Pdv >( aPropertyPointer );
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host::create_pdv(PDV aPdvType, moris::real aPdvVal)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.create_pdv() using pdv value with PDV type that doesn't exist on this host.");

            // create a pdv with pdv value
            mPdvList(mPdvTypeMap[aPdvType]) = std::make_shared< GEN_Pdv >( aPdvVal );
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Pdv_Host::is_active_type(PDV aPdvType)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.is_active_type() with type that doesn't exist on this host.");

            // return if active PDV
            return mActivePdvs(mPdvTypeMap[aPdvType]);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Pdv_Host::mark_pdv_as_inactive(PDV aPdvType)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                    "Tried to call Pdv_Host.mark_pdv_as_inactive() with type that doesn't exist on this host.");

            // Set PDV to be inactive
            mActivePdvs(mPdvTypeMap[aPdvType]) = false;
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Pdv_Host::get_global_index_for_pdv_type(PDV aPdvType)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                         "Tried to call Pdv_Host.get_global_index_for_pdv_type() with type that doesn't exist on this host.");

            // Return id from map
            return mGlobalPdvIndices(mPdvTypeMap[aPdvType]);
        }

        //--------------------------------------------------------------------------------------------------------------

        real Pdv_Host::get_pdv_value(PDV aPdvType)
        {
            // Check PDV type
            MORIS_ASSERT(mPdvTypeMap.key_exists(aPdvType),
                         "Tried to call Pdv_Host.get_pdv_value() with type that doesn't exist on this host.");

            // Return value
            return mPdvList(mPdvTypeMap[aPdvType])->get_val()(0);
        }

        //--------------------------------------------------------------------------------------------------------------
        
    } // end ge namepsace
} // end moris namespace

