/*
 * cl_FEM_Visualization_Set.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_FEM_VISUALIZATION_SET_HPP_
#define SRC_FEM_CL_FEM_VISUALIZATION_SET_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Visualization_Set.hpp"
#include "cl_FEM_Visualization_Object.hpp"
#include "cl_MTK_Set.hpp"

namespace moris
{
    namespace fem
    {
        class Vis_Set : public vis::Visualization_Set
        {
        private:

            moris::mtk::Set * mMeshSet = nullptr;

            Matrix< IndexMat > mElementIndexMap;
            uint mIndexOffset;

        public:
            Vis_Set( MSI::Equation_Set * aEquationSet,
                     moris::mtk::Set  * mMeshSet ) : Visualization_Set( aEquationSet ),
                                                     mMeshSet( mMeshSet )
            {
                this->create_visualization_objects();

//                this->create_index_map();
            };

//-----------------------------------------------------------------------------------------------------------

            ~Vis_Set()
            {

            };

//-----------------------------------------------------------------------------------------------------------

            Matrix< DDRMat > calculate_elemental_values()
            {
                MORIS_ERROR( mVisualizationObjects.size() != 0, "calculate_elemental_values(), no visualization objects. Did you call create_visualization_objects()");

                std::cout<<mElementIndexMap.numel()<<std::endl;
                std::cout<<"------------------------------------------------1-------"<<std::endl;

                Matrix< DDRMat > tMat( mElementIndexMap.numel(), 1, -1.0);
                for( auto tViosObj : mVisualizationObjects )
                {
                	Matrix< DDRMat > tVals = tViosObj->get_element_values();
                	Matrix< IndexMat > tEleInd = tViosObj->get_element_inds();

                	print( tVals, "tVals");
                	print( tEleInd, "tEleInd");
                	print( mElementIndexMap, "mElementIndexMap");
                    for( uint Ik = 0; Ik < tVals.numel(); Ik++ )
                    {
                        tMat( mElementIndexMap( tEleInd( Ik ) - mIndexOffset ) ) = tVals( Ik );
                    }
                }

                return tMat;
            }

//-----------------------------------------------------------------------------------------------------------

            void create_index_map()
            {
            	bool tBool = false;
                 mIndexOffset = mMeshSet->get_cell_inds_on_block( tBool ).min();

                 Matrix< IndexMat > tIndices = mMeshSet->get_cell_inds_on_block( tBool );

                 print(tIndices,"tIndices");

                 mElementIndexMap.resize( mMeshSet->get_cell_inds_on_block( tBool ).numel(),1 );

                 std::cout<<mElementIndexMap.numel()<<std::endl;
                 std::cout<<"------------------------------------------------1-------"<<std::endl;

                for( uint Ik = 0; Ik < mMeshSet->get_cell_inds_on_block( tBool ).numel(); Ik++ )
                {
                    std::cout<<tIndices( Ik ) - mIndexOffset<<std::endl;
                     mElementIndexMap( tIndices( Ik ) - mIndexOffset ) = Ik;
                }
            }

//-----------------------------------------------------------------------------------------------------------

            void create_visualization_objects()
            {
                uint tNumVisObjects = mEquationSet->get_num_equation_objects();

                mVisualizationObjects.resize( tNumVisObjects );

                for( uint Ik = 0; Ik < tNumVisObjects; Ik++ )
                {
                    mVisualizationObjects( Ik ) = new fem::Vis_Object( mEquationSet->get_equation_object_list()( Ik ),
                                                                       mMeshSet->get_clusters_by_index( Ik ) );
                }
            }

//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_VISUALIZATION_SET_HPP_ */
