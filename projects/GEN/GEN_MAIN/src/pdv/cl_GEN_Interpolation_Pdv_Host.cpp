#include "cl_GEN_Interpolation_Pdv_Host.hpp"
#include "cl_GEN_Pdv_Value.hpp"
#include "cl_GEN_Pdv_Property.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"

namespace moris
{
    namespace ge
    {
        
        //--------------------------------------------------------------------------------------------------------------
        
        Interpolation_Pdv_Host::Interpolation_Pdv_Host(
                Pdv_Host_Manager * aPdvHostManager,
                const moris_index & aNodeIndex,
                const moris_id & aNodeId,
                const moris_index & aNodeOwner,
                const Matrix<DDRMat> & aCoordinates,
                const Cell<PDV_Type> & aPDVTypes)
        : mPdvHostManager(aPdvHostManager),
          mNodeIndex(aNodeIndex),
          mNodeId(aNodeId),
          mNodeOwner(aNodeOwner),
          mCoordinates(aCoordinates),
          mPDVs(mPdvHostManager->get_max_num_pdvs(), nullptr)
        {

        }

        //--------------------------------------------------------------------------------------------------------------
        
        Interpolation_Pdv_Host::~Interpolation_Pdv_Host()
        {
        }

        //--------------------------------------------------------------------------------------------------------------
        
        uint Interpolation_Pdv_Host::get_num_pdvs()
        {
            return mPDVs.size();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::create_pdv(PDV_Type aPDVType, real aPdvVal)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            if( mPDVs(tPDVIndex) == nullptr )
            {
                // Create a pdv with pdv value
                mPDVs(tPDVIndex) = std::make_shared< Pdv_Value >(aPdvVal);
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::set_pdv_id(
                PDV_Type aPDVType,
                const moris_id aCounterId )
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            uint tPDVIndex = tPDVTypeMap( static_cast<uint>(aPDVType) );

            mPDVs( tPDVIndex )->set_id( aCounterId );
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Interpolation_Pdv_Host::get_pdv_id(PDV_Type aPDVType)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();
            uint tPDVIndex = tPDVTypeMap( static_cast<uint>(aPDVType) );

            return this->get_pdv_id(tPDVIndex);
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Interpolation_Pdv_Host::get_pdv_id(uint aPDVIndex)
        {
            return mPDVs(aPDVIndex)->get_id();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::create_pdv(PDV_Type aPDVType, std::shared_ptr<Property> aPropertyPointer)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // Create a pdv with property pointer
            mPDVs(tPDVIndex) = std::make_shared< Pdv_Property >(aPropertyPointer);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Interpolation_Pdv_Host::is_active_type(PDV_Type aPDVType)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // return if active PDV
            return mPDVs(tPDVIndex)->mIsActive;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::set_global_index_for_pdv_type(PDV_Type aPDVType, moris_id aId)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            mPDVs(tPDVIndex)->set_id( aId );
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Interpolation_Pdv_Host::get_global_index_for_pdv_type(PDV_Type aPDVType)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // Return id from map
            return mPDVs(tPDVIndex)->get_id();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Interpolation_Pdv_Host::get_all_global_indices()
        {
            // Initialize global IDs
            Matrix<DDSMat> tGlobPdvId(mPDVs.size(), 1);

            // Loop over PDVs and get IDs
            uint tCounter = 0;
            for (auto tPdv : mPDVs)
            {
                if (tPdv != nullptr)
                {
                    tGlobPdvId(tCounter++) = tPdv->get_id();
                }
            }

            // Resize IDs
            tGlobPdvId.resize(tCounter, 1);

            // Return IDs
            return tGlobPdvId;
        }

        //--------------------------------------------------------------------------------------------------------------

        real Interpolation_Pdv_Host::get_pdv_value(PDV_Type aPDVType)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // Return value
            return mPDVs(tPDVIndex)->get_value(mNodeIndex, mCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Interpolation_Pdv_Host::get_pdv_exists(PDV_Type aPDVType)
        {
            const Matrix< DDSMat > & tPDVTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPDVIndex = tPDVTypeMap( static_cast<sint>(aPDVType) );

            // Check PDV type
            MORIS_ASSERT(tPDVIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            bool tBool = false;

            if( mPDVs(tPDVIndex) != nullptr )
            {
                tBool = true;
            }

            return tBool;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Interpolation_Pdv_Host::get_sensitivities(uint aPDVIndex)
        {
            // If PDV exists and is active, ask it for sensitivities. Otherwise, return nothing
            return ((mPDVs(aPDVIndex) and mPDVs(aPDVIndex)->mIsActive)
                    ? mPDVs(aPDVIndex)->get_sensitivities(mNodeIndex, mCoordinates) : Matrix<DDRMat>(0, 0));
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Interpolation_Pdv_Host::get_determining_adv_ids(uint aPDVIndex)
        {
            return mPDVs(aPDVIndex)->get_determining_adv_ids(mNodeIndex, mCoordinates);
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}
