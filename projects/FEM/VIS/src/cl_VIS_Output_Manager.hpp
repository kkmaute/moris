/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Output_Manager.hpp
 *
 */

#ifndef SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_
#define SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_

#include "cl_Vector.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Output_Enums.hpp"

#include "cl_MSI_Equation_Set.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_Parameter_List.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh_Manager;
    }

    namespace MSI
    {
        class Equation_Model;
    }

    namespace mdl
    {
        class Model;
    }

    namespace vis
    {
        enum class File_Status
        {
            CLOSED,
            OPEN,
            CLOSINGDELAYED,
            END_ENUM
        };

        struct Output_Data
        {
            //! Mesh Index for sanity checks
            sint mMeshIndex = -1;

            //! Frequency with which output file is save during transient simulation
            sint mSaveFrequency = MORIS_SINT_MAX;

            //! Counter of writing fields to mesh
            sint mFieldWriteCounter = 0;

            //! Time offset for writing sequence of optimization steps
            real mTimeOffset = 0.0;

            //! Mesh Type
            enum VIS_Mesh_Type mMeshType;

            //! flag whether IQIs exist that are based on sensitivity analysis
            bool mSensitivityBasedIQIs;

            //! Output Path
            std::string mOutputPath;

            //! Mesh Name
            std::string mMeshName;

            //! Output Path for temporary file
            std::string mTempPath;

            //! Temporary Name
            std::string mTempName;

            //! Mesh Name
            std::string mMeshPath;

            //! Set names which shall be part of this mesh
            Vector< std::string > mSetNames;

            //! Field names which shall be used for outputs
            Vector< std::string > mFieldNames;

            //! Field types
            Vector< enum Field_Type > mFieldType;

            //! Quantity of interest names
            Vector< std::string > mQINames;

            //! Analysis type for which quantity of interest is evaluated
            Vector< enum Analysis_Type > mAnalysisType;
        };

        //-----------------------------------------------------------------------------------------------------------

        class Output_Manager
        {
          private:
            Vector< vis::Output_Data > mOutputData;

            Vector< mtk::Mesh * > mVisMesh;

            //! file closing status
            Vector< vis::File_Status > mVisMeshOutputFileStatus;

            bool mOnlyPrimary = false;

            Vector< moris::mtk::Writer_Exodus * > mWriter;

            std::shared_ptr< mtk::Mesh_Manager > mMTKMesh = nullptr;

            moris::uint mMTKMeshPairIndex = MORIS_UINT_MAX;

            moris::real mTimeShift = 0.0;

          protected:

          public:
            Output_Manager() {};

            //-----------------------------------------------------------------------------------------------------------

            Output_Manager( const moris::Parameter_List &aParameterlist )
            {
                this->set_outputs( aParameterlist );
            };

            //-----------------------------------------------------------------------------------------------------------

            Output_Manager( Vector< moris::Parameter_List > aParameterList )
            {
                for ( uint iVisParameter = 0; iVisParameter < aParameterList.size(); iVisParameter++ )
                {
                    this->set_outputs( aParameterList( iVisParameter ) );
                }
            }

            //-----------------------------------------------------------------------------------------------------------

            ~Output_Manager()
            {
                // close files
                for ( uint i = 0; i < mVisMeshOutputFileStatus.size(); ++i )
                {
                    if ( mVisMeshOutputFileStatus( i ) != File_Status::CLOSED )
                    {
                        this->end_writing( i );
                    }
                }

                // delete vis meshes
                for ( auto tMesh : mVisMesh )
                {
                    delete tMesh;
                }
            };

            //-----------------------------------------------------------------------------------------------------------

            void
            delete_pointers( const uint aVisMeshIndex )
            {
                delete mVisMesh( aVisMeshIndex );

                mVisMesh( aVisMeshIndex ) = nullptr;

                delete mWriter( aVisMeshIndex );

                mWriter( aVisMeshIndex ) = nullptr;
            }

            //-----------------------------------------------------------------------------------------------------------

            void
            end_writing( const uint aVisMeshIndex )
            {
                // check if file is already closed
                MORIS_ERROR( mVisMeshOutputFileStatus( aVisMeshIndex ) != File_Status::CLOSED,
                        "Output_Manager::end_writing - closing of already closed file not permitted" );

                // only close output file if mesh is not empty
                if ( mWriter( aVisMeshIndex ) != nullptr )
                {
                    mWriter( aVisMeshIndex )->close_file();
                }

                this->delete_pointers( aVisMeshIndex );

                // set file closing status to closed
                mVisMeshOutputFileStatus( aVisMeshIndex ) = File_Status::CLOSED;
            }

            //-----------------------------------------------------------------------------------------------------------

            bool
            check_for_closing_file( uint aVisMeshIndex, bool aCloseFile, bool aIsFowardAnalysis )
            {
                // check if requested action is to not close file
                if ( !aCloseFile )
                {
                    return false;
                }

                // check if file is already closed
                if ( mVisMeshOutputFileStatus( aVisMeshIndex ) == File_Status::CLOSED )
                {
                    MORIS_ERROR( false,
                            "Output_Manager::check_for_closing_file - cannot close file whose status is closed" );

                    return false;
                }

                // if in forward analysis
                if ( aIsFowardAnalysis )
                {
                    // file status is closing delayed, close file (can happen in sweep)
                    if ( mVisMeshOutputFileStatus( aVisMeshIndex ) == File_Status::CLOSINGDELAYED )
                    {
                        MORIS_ERROR( false,
                                "Output_Manager::check_for_closing_file - file status is closing delayed in forward analysis" );
                    }

                    // sensitivity based IQIs exist, delay closing
                    if ( mOutputData( aVisMeshIndex ).mSensitivityBasedIQIs )
                    {
                        mVisMeshOutputFileStatus( aVisMeshIndex ) = File_Status::CLOSINGDELAYED;
                        return false;
                    }
                }
                // if in sensitivity analysis, close file
                else
                {
                    // file status is closing delayed, keep file open and file status to open
                    if ( mVisMeshOutputFileStatus( aVisMeshIndex ) == File_Status::CLOSINGDELAYED )
                    {
                        mVisMeshOutputFileStatus( aVisMeshIndex ) = File_Status::OPEN;

                        return true;
                    }
                    else
                    {
                        // file status is open, close file
                        MORIS_ERROR( false,
                                "Output_Manager::check_for_closing_file - in sensitivity analysis file status should be closing delayed" );

                        return true;
                    }
                }

                // for all other cases close file
                return true;
            }

            //-----------------------------------------------------------------------------------------------------------

            void set_outputs(
                    const uint                          aOutputIndex,
                    const enum VIS_Mesh_Type            aMeshType,
                    const std::string                  &aMeshPath,
                    const std::string                  &aMeshName,
                    const std::string                  &aTempPath,
                    const std::string                  &aTempName,
                    const Vector< std::string >        &aBlockNames,
                    const Vector< std::string >        &aFieldNames,
                    const Vector< enum Field_Type >    &aFieldType,
                    const Vector< std::string >        &aQINames,
                    const Vector< enum Analysis_Type > &aAnalysisType,
                    const uint                          aSaveFrequency = 1,
                    const real                          aTimeOffset    = 0.0 );

            //---------------------------------------------------------------------------------------------------------------------------

            void set_outputs( const moris::Parameter_List &aParameterlist );

            //---------------------------------------------------------------------------------------------------------------------------

            void setup_vis_mesh_for_output(
                    const uint                                    aVisMeshIndex,
                    const std::shared_ptr< mtk::Mesh_Manager >   &aMesh,
                    const uint                                    aMeshPairIndex,
                    const std::shared_ptr< MSI::Equation_Model > &aEquationModel );

            //---------------------------------------------------------------------------------------------------------------------------

            void create_visualization_mesh(
                    const uint                                  aVisMeshIndex,
                    const std::shared_ptr< mtk::Mesh_Manager > &aMesh,
                    const uint                                  aMeshPairIndex );

            //-----------------------------------------------------------------------------------------------------------

            void set_visualization_sets(
                    const uint                                    aVisMeshIndex,
                    const std::shared_ptr< MSI::Equation_Model > &aEquationModel );

            //-----------------------------------------------------------------------------------------------------------

            void write_mesh( const uint aVisMeshIndex );

            //-----------------------------------------------------------------------------------------------------------

            void write_mesh_indices( const uint aVisMeshIndex );

            //-----------------------------------------------------------------------------------------------------------

            void add_nodal_fields( const uint aVisMeshIndex );

            //-----------------------------------------------------------------------------------------------------------

            void add_elemental_fields( const uint aVisMeshIndex );

            //-----------------------------------------------------------------------------------------------------------

            void add_faceted_fields( const uint aVisMeshIndex );

            //-----------------------------------------------------------------------------------------------------------

            void add_global_fields( const uint aVisMeshIndex );

            //-----------------------------------------------------------------------------------------------------------

            void write_field(
                    const uint                                    aVisMeshIndex,
                    const real                                    aTime,
                    const std::shared_ptr< MSI::Equation_Model > &aEquationModel,
                    const bool                                    aIsFowardAnalysis );

            //-----------------------------------------------------------------------------------------------------------

            void
            get_IQI_and_field_names(
                    const uint                       aVisMeshIndex,
                    Vector< Vector< std::string > > &aIQINames,
                    Vector< Vector< std::string > > &aFieldNames,
                    Vector< uint >                  &aNumIQIsForFieldType,
                    const bool                       aIsFowardAnalysis );

            //-----------------------------------------------------------------------------------------------------------

            void compute_fields_for_set(
                    const uint                             aVisMeshIndex,
                    MSI::Equation_Set                     *aFemSet,
                    Vector< Vector< std::string > > const &aIQINames,
                    Vector< Vector< std::string > > const &aFieldNames,
                    Matrix< DDRMat >                      *aGlobalFieldValues,
                    Matrix< DDRMat >                      *aNodalFieldValues );

            //-----------------------------------------------------------------------------------------------------------

            void
            compute_and_write_elemental_fields_on_set(
                    const uint                   aVisMeshIndex,
                    MSI::Equation_Set           *aFemSet,
                    const Field_Type             aFieldType,
                    Vector< std::string > const &aIQINamesForType,
                    Vector< std::string > const &aFieldNamesForType );

            //-----------------------------------------------------------------------------------------------------------

            bool sensitivity_based_iqis_exist( uint aVisMeshIndex )
            {
                return mOutputData( aVisMeshIndex ).mSensitivityBasedIQIs;
            }

            //-----------------------------------------------------------------------------------------------------------
        };
    }    // namespace vis
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_ */
