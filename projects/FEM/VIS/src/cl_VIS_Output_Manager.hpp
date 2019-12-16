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

namespace moris
{
    namespace mtk
    {
        class Mesh_Manager;
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

            //! Mesh Name
            std::string                     mMeshName;

            //! Set names which shall be part of this mesh
            moris::Cell< std::string >      mSetNames;

            //! Set indices
            moris::Cell< moris_index >      mSetIndices;

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

            moris::Cell< MSI::Equation_Set * > mEquationSets;

            Cell< mtk::Mesh * > mVisMesh;

            Matrix< IndexMat > mListOfRequestedBlocks;

            bool mOnlyPrimary = false;

            moris::Cell< Writer_Exodus * >  mWriter;

        protected:

        public:
            Output_Manager()
            {
                mListOfRequestedBlocks = { { 0, 1, 2, 3 } };
                mOnlyPrimary = false ;

            };

//-----------------------------------------------------------------------------------------------------------

            ~Output_Manager()
            {
               //FIXME add delete
            };

//-----------------------------------------------------------------------------------------------------------

            void delete_pointers( const uint aVisMeshIndex )
            {
                delete ( mVisMesh( aVisMeshIndex ) );

                mVisMesh( aVisMeshIndex ) = nullptr;

                delete ( mWriter( aVisMeshIndex ) );

                mWriter( aVisMeshIndex ) = nullptr;
            }

//-----------------------------------------------------------------------------------------------------------

            void end_writing( const uint aVisMeshIndex )
            {
                mWriter( aVisMeshIndex )->close_file();

                this->delete_pointers( aVisMeshIndex );
            }

//-----------------------------------------------------------------------------------------------------------

            void set_outputs( const uint                              aOutputIndex,
                              const enum VIS_Mesh_Type                aMeshType,
                              const std::string                     & aMeshName,
                              const moris::Cell< std::string >      & aBlockNames,
                              const moris::Cell< moris_index >      & aBlockIndices,
                              const moris::Cell< std::string >      & aFieldNames,
                              const moris::Cell< enum Field_Type >  & aFieldType,
                              const moris::Cell< enum Output_Type > & aEnum );

//---------------------------------------------------------------------------------------------------------------------------

            void create_visualization_mesh( const uint                aVisMeshIndex,
                                                  mtk::Mesh_Manager * aMesh,
                                            const uint                aMeshPairIndex);

//-----------------------------------------------------------------------------------------------------------

            void set_visualization_sets( const uint         aVisMeshIndex,
                                               mdl::Model * aModel );

//-----------------------------------------------------------------------------------------------------------

            void write_mesh( const uint aVisMeshIndex,
                             const real tTime);

//-----------------------------------------------------------------------------------------------------------

            void write_mesh_indices( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void add_nodal_fields( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void add_elemetal_fields( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void add_global_fields( const uint aVisMeshIndex );

//-----------------------------------------------------------------------------------------------------------

            void write_field( const uint         aVisMeshIndex,
                                    mdl::Model * aModel );

//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_ */
