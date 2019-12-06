/*
 * cl_VIS_Output_Data.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_
#define SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Factory.hpp"
#include "cl_FEM_Visualization_Set.hpp"

#include "cl_MSI_Equation_Set.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

namespace moris
{
    namespace vis
    {
        class Output_Data
        {
        private:
            moris::Cell< MSI::Equation_Set * > & mEquationSets;

            mtk::Mesh_Manager         * mMesh;

            const uint                  mMeshPairIndex;

            Cell< mtk::Mesh * > mVisMesh;

            Cell< Cell< Visualization_Set * > > mVisualizationSets;

            Matrix< IndexMat > mListOfRequestedBlocks;

            bool mOnlyPrimary = false;

            Writer_Exodus *  mWriter =nullptr;

        protected:

        public:
            Output_Data(       moris::Cell< MSI::Equation_Set * > & aEquationSets,
                               mtk::Mesh_Manager                  * aMesh,
                         const uint                                 aMeshPairIndex ) : mEquationSets( aEquationSets ),
                                                                                       mMesh( aMesh ),
                                                                                       mMeshPairIndex( aMeshPairIndex )
            {
                mListOfRequestedBlocks = { { 0, 1, 2, 3 } };
                mOnlyPrimary = true;

                this->create_visualization_meshes();

                this->create_visualization_sets();

                mWriter = new Writer_Exodus(mVisMesh( 0 ));

            };

//-----------------------------------------------------------------------------------------------------------

            ~Output_Data()
            {
               //FIXME add delete
            };

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_meshes()
            {
                vis::Factory tVisFactory( mMesh, mMeshPairIndex );

                mVisMesh.resize( 1, nullptr );

                for( uint Ik = 0; Ik < 1; Ik++ )
                {
                    mVisMesh( Ik ) = tVisFactory.create_visualization_mesh();
                }
            }

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_sets()
            {
                mVisualizationSets.resize( 1 );

                for( uint Ik = 0; Ik < 1; Ik++ )
                {
                    uint tRequestedSets = mListOfRequestedBlocks.numel();
                    mVisualizationSets( Ik ).resize( tRequestedSets );

                    for( uint Ii = 0; Ii < tRequestedSets; Ii++ )
                    {
                        uint tBlockIndex = mListOfRequestedBlocks( Ii );
                        mVisualizationSets( Ik )( Ii ) = new fem::Vis_Set( mEquationSets( tBlockIndex ),
                                                                           mVisMesh( Ik )->get_block_by_index( tBlockIndex ) );
                    }
                }
            }

//-----------------------------------------------------------------------------------------------------------

            void write_mesh()
            {
                std::string tPrefix = std::getenv("MORISROOT");
                std::string tMeshFilePath = tPrefix + "build";
            	mWriter->write_mesh(tMeshFilePath, "Vis_Mesh_3.exo");

            	mWriter->close_file();
            }

//-----------------------------------------------------------------------------------------------------------

            void write_field()
            {
                std::string tPrefix = std::getenv("MORISROOT");
                std::string tMeshFileName = tPrefix + "build/Vis_Mesh_3.exo";
            	mWriter->open_file(tMeshFileName, false );

                moris::Cell<std::string> tElementalFieldNames(1);
                tElementalFieldNames(0) = "pressure";

                //-------------------------------------------------------------------------------------------

//                moris::Cell<const moris::mtk::Cell*> tElementsInBlock = mVisMesh( 0 )->get_block_set_cells("HMR_dummy_c_p0");
//
//                uint tNumElements = tElementsInBlock.size();
//                moris::Matrix<moris::DDRMat> tetField(tNumElements, 1, 4);
//
//
//                for(uint Ik = 0; Ik<tNumElements;Ik++)
//                {
//                	tetField( Ik ) = Ik;
//                }

                //-------------------------------------------------------------------------------------------
                for ( uint Ik = 0; Ik<4;Ik++ )
                {
                    moris::Matrix<moris::DDRMat> tetField = mVisualizationSets( 0 )( Ik )->calculate_elemental_values();

                    mWriter->set_elemental_fields(tElementalFieldNames);
                    mWriter->set_time(0.0);
                    mWriter->write_elemental_field( Ik, "pressure", tetField);
                }

//                mWriter->set_elemental_fields(tElementalFieldNames);
//                mWriter->set_time(0.0);
//                mWriter->write_elemental_field( 0, "pressure", tetField);

                mWriter->close_file();
            }


//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_OUTPUT_DATA_HPP_ */
