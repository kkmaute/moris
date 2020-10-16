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
                const Cell<PDV_Type> & aPdvTypes)
        : mPdvHostManager(aPdvHostManager),
          mNodeIndex(aNodeIndex),
          mNodeId(aNodeId),
          mNodeOwner(aNodeOwner),
          mCoordinates(aCoordinates),
          mPdvs(mPdvHostManager->get_max_num_pdvs(), nullptr)
        {

        }

        //--------------------------------------------------------------------------------------------------------------
        
        Interpolation_Pdv_Host::~Interpolation_Pdv_Host()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::create_pdv(PDV_Type aPdvType, real aPdvVal)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            if( mPdvs(tPdvHostIndex) == nullptr )
            {
                // Create a pdv with pdv value
                mPdvs(tPdvHostIndex) = std::make_shared< Pdv_Value >(aPdvVal);
            }
            else
            {
                mPdvs(tPdvHostIndex)->set_value( aPdvVal );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::set_pdv_id(
                enum PDV_Type aPdvType,
                const moris_id aCounterId )
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            uint tPdvHostIndex = tPdvTypeMap( static_cast<uint>(aPdvType) );

            mPdvs( tPdvHostIndex )->set_id( aCounterId );
        }

        //--------------------------------------------------------------------------------------------------------------

        moris_id Interpolation_Pdv_Host::get_pdv_id( enum PDV_Type aPdvType)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            uint tPdvHostIndex = tPdvTypeMap( static_cast<uint>(aPdvType) );

            return mPdvs( tPdvHostIndex )->get_id();
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::create_pdv(PDV_Type aPdvType, std::shared_ptr<Property> aPropertyPointer)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // Create a pdv with property pointer
            mPdvs(tPdvHostIndex) = std::make_shared< Pdv_Property >(aPropertyPointer);
        }

        //--------------------------------------------------------------------------------------------------------------

        bool Interpolation_Pdv_Host::is_active_type(PDV_Type aPdvType)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // return if active PDV
            return mPdvs(tPdvHostIndex)->mIsActive;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::set_global_index_for_pdv_type(PDV_Type aPdvType, moris_id aId)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            mPdvs(tPdvHostIndex)->set_id( aId );
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Interpolation_Pdv_Host::get_global_index_for_pdv_type(PDV_Type aPdvType)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // Return id from map
            return mPdvs(tPdvHostIndex)->get_id();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDSMat> Interpolation_Pdv_Host::get_all_global_indices()
        {
            // Initialize global IDs
            Matrix<DDSMat> tGlobPdvId(mPdvs.size(), 1);

            // Loop over PDVs and get IDs
            uint tCounter = 0;
            for (auto tPdv : mPdvs)
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

        real Interpolation_Pdv_Host::get_pdv_value(PDV_Type aPdvType)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            // Return value
            return mPdvs(tPdvHostIndex)->get_value(mNodeIndex, mCoordinates);
        }

        bool Interpolation_Pdv_Host::get_pdv_exists(PDV_Type aPdvType)
        {
            const Matrix< DDSMat > & tPdvTypeMap = mPdvHostManager->get_pdv_type_map();

            sint tPdvHostIndex = tPdvTypeMap( static_cast<sint>(aPdvType) );

            // Check PDV type
            MORIS_ASSERT(tPdvHostIndex != -1,
                         "Tried to call Pdv_Host.create_pdv() using pdv value with a PDV type that doesn't exist on this host.");

            bool tBool = false;

            if( mPdvs(tPdvHostIndex) != nullptr )
            {
                tBool = true;
            }

            return tBool;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Interpolation_Pdv_Host::get_all_sensitivities(Matrix<DDRMat>& aSensitivities)
        {
            aSensitivities.resize(0, aSensitivities.n_cols());
            if (mPdvs.size() > 0)
            {
                uint tActivePdvIndex = 0;
                Matrix<DDRMat> tSensitivities(0, 0);
                for (uint tPdvIndex = 0; tPdvIndex < mPdvs.size(); tPdvIndex++)
                {
                    if( mPdvs(tPdvIndex) != nullptr)
                    {
                        if (mPdvs(tPdvIndex)->mIsActive)
                        {
                            mPdvs(tPdvIndex)->get_sensitivity(mNodeIndex, mCoordinates, tSensitivities);
                            aSensitivities.resize(tActivePdvIndex + 1, tSensitivities.length());
                            aSensitivities.set_row(tActivePdvIndex++, tSensitivities);
                        }
                    }
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        
    }
}
