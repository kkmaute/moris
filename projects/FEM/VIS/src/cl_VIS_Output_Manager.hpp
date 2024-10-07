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
            Vector< std::string > mSetNames;

            //! Field names which shall be used for outputs
            Vector< std::string > mFieldNames;

            //! Field types
            Vector< enum Field_Type > mFieldType;

            //! Quantity of interest names
            Vector< std::string > mQINames;
        };

        //-----------------------------------------------------------------------------------------------------------

        class Output_Manager
        {
          private:
            Vector< vis::Output_Data > mOutputData;

            Vector< mtk::Mesh * > mVisMesh;
            Vector< bool >        mVisMeshCreatedAndOpen;

            bool mOnlyPrimary = false;

            Vector< moris::mtk::Writer_Exodus * > mWriter;

            std::shared_ptr< mtk::Mesh_Manager > mMTKMesh = nullptr;

            moris::uint mMTKMeshPairIndex = MORIS_UINT_MAX;

            moris::real mTimeShift = 0.0;

          protected:

          public:
            Output_Manager(){};

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
                    const Vector< std::string >     &aBlockNames,
                    const Vector< std::string >     &aFieldNames,
                    const Vector< enum Field_Type > &aFieldType,
                    const Vector< std::string >     &aQINames,
                    const uint                            aSaveFrequency = 1,
                    const real                            aTimeOffset    = 0.0 );

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
                    const std::shared_ptr< MSI::Equation_Model > &aEquationModel );

            //-----------------------------------------------------------------------------------------------------------

            void
            get_IQI_and_field_names(
                    const uint                   aVisMeshIndex,
                    Vector< Vector< std::string > > &aIQINames,
                    Vector< Vector< std::string > > &aFieldNames,
                    Vector< uint >                &aNumIQIsForFieldType );

            //-----------------------------------------------------------------------------------------------------------

            void compute_fields_for_set(
                    const uint                         aVisMeshIndex,
                    MSI::Equation_Set                 *aFemSet,
                    Vector< Vector< std::string > > const &aIQINames,
                    Vector< Vector< std::string > > const &aFieldNames,
                    Matrix< DDRMat >                  *aGlobalFieldValues,
                    Matrix< DDRMat >                  *aNodalFieldValues );

            //-----------------------------------------------------------------------------------------------------------

            void
            compute_and_write_elemental_fields_on_set(
                    const uint                 aVisMeshIndex,
                    MSI::Equation_Set         *aFemSet,
                    const Field_Type           aFieldType,
                    Vector< std::string > const &aIQINamesForType,
                    Vector< std::string > const &aFieldNamesForType );

            //-----------------------------------------------------------------------------------------------------------
        };
    }    // namespace vis
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_ */
