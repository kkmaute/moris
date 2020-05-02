/*
 * cl_VIS_Output_Manager.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
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
            sint                            mMeshIndex = -1;

            //! Mesh Type
            enum VIS_Mesh_Type              mMeshType;

            //! Output Path
            std::string                     mOutputPath;

            //! Mesh Name
            std::string                     mMeshName;

            //! Mesh Name
            std::string                     mMeshPath;

            //! Set names which shall be part of this mesh
            moris::Cell< std::string >      mSetNames;

            //! Field names which shall be used for outputs
            moris::Cell< std::string >      mFieldNames;

            //! Field types
            moris::Cell< enum Field_Type >  mFieldType;

            //! Output enums
            moris::Cell< enum Output_Type > mOutputType;
        };

//-----------------------------------------------------------------------------------------------------------

        class Output_Manager
        {
        private:

            moris::Cell< vis::Output_Data > mOutputData;

            moris::Cell< mtk::Mesh * > mVisMesh;
            moris::Cell< bool >        mVisMeshCreatedAndOpen;

            bool mOnlyPrimary = false;

            moris::Cell< moris::mtk::Writer_Exodus * >  mWriter;

            mtk::Mesh_Manager *             mMTKMesh = nullptr;

            moris::uint                     mMTKMeshPairIndex;

            moris::real mTimeStamp;

        protected:

        public:
            Output_Manager(){};

//-----------------------------------------------------------------------------------------------------------

            Output_Manager( moris::ParameterList aParamterelist )
            {
                this->set_outputs( aParamterelist );
            };

//-----------------------------------------------------------------------------------------------------------

            ~Output_Manager()
            {
                for( auto tMesh : mVisMesh )
                {
                    delete tMesh;
                }
            };

//-----------------------------------------------------------------------------------------------------------

            void delete_pointers( const uint aVisMeshIndex )
            {
                delete ( mVisMesh( aVisMeshIndex ) );

                mVisMesh( aVisMeshIndex ) = nullptr;

                delete ( mWriter( aVisMeshIndex ) );

                mWriter( aVisMeshIndex ) = nullptr;

                mVisMeshCreatedAndOpen( aVisMeshIndex ) = false;
            }

//-----------------------------------------------------------------------------------------------------------

            void end_writing( const uint aVisMeshIndex )
            {
                mWriter( aVisMeshIndex )->close_file();

                this->delete_pointers( aVisMeshIndex );

                MORIS_LOG( " Finished writing output. ");
            }

//-----------------------------------------------------------------------------------------------------------

            void set_outputs( const uint                              aOutputIndex,
                              const enum VIS_Mesh_Type                aMeshType,
                              const std::string                     & aMeshPath,
                              const std::string                     & aMeshName,
                              const moris::Cell< std::string >      & aBlockNames,
                              const moris::Cell< std::string >      & aFieldNames,
                              const moris::Cell< enum Field_Type >  & aFieldType,
                              const moris::Cell< enum Output_Type > & aEnum );

//---------------------------------------------------------------------------------------------------------------------------

            void set_outputs( moris::ParameterList aParamterelist );

//---------------------------------------------------------------------------------------------------------------------------

            void create_visualization_mesh( const uint                aVisMeshIndex,
                                                  mtk::Mesh_Manager * aMesh,
                                            const uint                aMeshPairIndex);

//-----------------------------------------------------------------------------------------------------------

            void set_visualization_sets( const uint                                   aVisMeshIndex,
                                               std::shared_ptr< MSI::Equation_Model > aEquationModel );

//-----------------------------------------------------------------------------------------------------------

            void write_mesh( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void write_mesh_indices( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void add_nodal_fields( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void add_elemetal_fields( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void add_global_fields( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void write_field( const uint                                   aVisMeshIndex,
                              const real                                   aTime,
                                    std::shared_ptr< MSI::Equation_Model > aEquationModel );

//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_ */
