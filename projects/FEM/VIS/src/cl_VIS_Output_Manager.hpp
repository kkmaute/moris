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

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Output_Enums.hpp"

#include "cl_MSI_Equation_Set.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_Param_List.hpp"

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
            moris::Cell< std::string > mSetNames;

            //! Field names which shall be used for outputs
            moris::Cell< std::string > mFieldNames;

            //! Field types
            moris::Cell< enum Field_Type > mFieldType;

            //! Quantity of interest names
            moris::Cell< std::string > mQINames;
        };

        //-----------------------------------------------------------------------------------------------------------

        class Output_Manager
        {
          private:
            moris::Cell< vis::Output_Data > mOutputData;

            moris::Cell< mtk::Mesh * > mVisMesh;
            moris::Cell< bool >        mVisMeshCreatedAndOpen;

            bool mOnlyPrimary = false;

            moris::Cell< moris::mtk::Writer_Exodus * > mWriter;

            std::shared_ptr< mtk::Mesh_Manager > mMTKMesh = nullptr;

            moris::uint mMTKMeshPairIndex = MORIS_UINT_MAX;

            moris::real mTimeShift = 0.0;

          protected:

          public:
            Output_Manager(){};

            //-----------------------------------------------------------------------------------------------------------

            Output_Manager( moris::ParameterList aParameterlist )
            {
                this->set_outputs( aParameterlist );
            };

            //-----------------------------------------------------------------------------------------------------------

            Output_Manager( moris::Cell< moris::ParameterList > aParameterList )
            {
                for ( uint iVisParameter = 0; iVisParameter < aParameterList.size(); iVisParameter++ )
                {
                    this->set_outputs( aParameterList( iVisParameter ) );
                }
            }

            //-----------------------------------------------------------------------------------------------------------

            ~Output_Manager()
            {
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

                mVisMeshCreatedAndOpen( aVisMeshIndex ) = false;
            }

            //-----------------------------------------------------------------------------------------------------------

            void
            end_writing( const uint aVisMeshIndex )
            {
                // only close output file if mesh is not empty
                if ( mWriter( aVisMeshIndex ) != nullptr )
                {
                    mWriter( aVisMeshIndex )->close_file();
                }

                this->delete_pointers( aVisMeshIndex );
            }

            //-----------------------------------------------------------------------------------------------------------

            void set_outputs(
                    const uint                            aOutputIndex,
                    const enum VIS_Mesh_Type              aMeshType,
                    const std::string                    &aMeshPath,
                    const std::string                    &aMeshName,
                    const std::string                    &aTempPath,
                    const std::string                    &aTempName,
                    const moris::Cell< std::string >     &aBlockNames,
                    const moris::Cell< std::string >     &aFieldNames,
                    const moris::Cell< enum Field_Type > &aFieldType,
                    const moris::Cell< std::string >     &aQINames,
                    const uint                            aSaveFrequency = 1,
                    const real                            aTimeOffset    = 0.0 );

            //---------------------------------------------------------------------------------------------------------------------------

            void set_outputs( moris::ParameterList aParameterlist );

            //---------------------------------------------------------------------------------------------------------------------------

            void setup_vis_mesh_for_output(
                    const uint                             aVisMeshIndex,
                    std::shared_ptr< mtk::Mesh_Manager >   aMesh,
                    const uint                             aMeshPairIndex,
                    std::shared_ptr< MSI::Equation_Model > aEquationModel );

            //---------------------------------------------------------------------------------------------------------------------------

            void create_visualization_mesh(
                    const uint                           aVisMeshIndex,
                    std::shared_ptr< mtk::Mesh_Manager > aMesh,
                    const uint                           aMeshPairIndex );

            //-----------------------------------------------------------------------------------------------------------

            void set_visualization_sets(
                    const uint                             aVisMeshIndex,
                    std::shared_ptr< MSI::Equation_Model > aEquationModel );

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
                    const uint                             aVisMeshIndex,
                    const real                             aTime,
                    std::shared_ptr< MSI::Equation_Model > aEquationModel );

            //-----------------------------------------------------------------------------------------------------------

            void
            get_IQI_and_field_names(
                    const uint                   aVisMeshIndex,
                    Cell< Cell< std::string > > &aIQINames,
                    Cell< Cell< std::string > > &aFieldNames,
                    Cell< uint >                &aNumIQIsForFieldType );

            //-----------------------------------------------------------------------------------------------------------

            void compute_fields_for_set(
                    const uint                         aVisMeshIndex,
                    MSI::Equation_Set                 *aFemSet,
                    Cell< Cell< std::string > > const &aIQINames,
                    Cell< Cell< std::string > > const &aFieldNames,
                    Matrix< DDRMat >                  *aGlobalFieldValues,
                    Matrix< DDRMat >                  *aNodalFieldValues );

            //-----------------------------------------------------------------------------------------------------------

            void
            compute_and_write_elemental_fields_on_set(
                    const uint                 aVisMeshIndex,
                    MSI::Equation_Set         *aFemSet,
                    const Field_Type           aFieldType,
                    Cell< std::string > const &aIQINamesForType,
                    Cell< std::string > const &aFieldNamesForType );

            //-----------------------------------------------------------------------------------------------------------
        };
    }    // namespace vis
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_ */
