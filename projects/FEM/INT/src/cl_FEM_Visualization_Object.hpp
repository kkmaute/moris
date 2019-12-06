/*
 * cl_FEM_Visualization_Object.hpp
 *
 *  Created on: Dez 02, 2019
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_FEM_VISUALIZATION_OBJECT_HPP_
#define SRC_FEM_CL_FEM_VISUALIZATION_OBJECT_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_VIS_Visualization_Object.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris
{
    namespace fem
    {
        class Vis_Object : public vis::Visualization_Object
        {
        private:
            const mtk::Cluster * mMeshCluster = nullptr;


        public:
            Vis_Object(       MSI::Equation_Object * aEquationObject,
                        const mtk::Cluster         * aMeshCluster) : Visualization_Object( aEquationObject ),
                                                                     mMeshCluster ( aMeshCluster )
            {

            };

//-----------------------------------------------------------------------------------------------------------

            ~Vis_Object()
            {

            };

//-----------------------------------------------------------------------------------------------------------

            Matrix< DDRMat > get_element_values()
		    {
            	uint tnumPrimCells = mMeshCluster->get_num_primary_cells();
            	uint tnumVoidCells = mMeshCluster->get_num_void_cells();

            	Matrix< DDRMat > tVals( tnumPrimCells + tnumVoidCells, 1, -1.0 );

            	uint tCounter = 0;

                for( uint Ik = 0; Ik < tnumPrimCells; Ik++ )
                {
                	tVals( tCounter++) = mMeshCluster->get_primary_cells_in_cluster()( Ik )->get_id();
                }

                for( uint Ik = 0; Ik < tnumVoidCells; Ik++ )
                {
                	tVals( tCounter++) = mMeshCluster->get_void_cells_in_cluster()( Ik )->get_id();
                }
            	return tVals;
		    }

//-----------------------------------------------------------------------------------------------------------

            Matrix< IndexMat > get_element_inds()
		    {
            	uint tnumPrimCells = mMeshCluster->get_num_primary_cells();
            	uint tnumVoidCells = mMeshCluster->get_num_void_cells();

            	Matrix< IndexMat > tVals( tnumPrimCells + tnumVoidCells, 1  );

            	uint tCounter = 0;

                for( uint Ik = 0; Ik < tnumPrimCells; Ik++ )
                {
                	tVals( tCounter++) = mMeshCluster->get_primary_cells_in_cluster()( Ik )->get_index();
                }

                for( uint Ik = 0; Ik < tnumVoidCells; Ik++ )
                {
                	tVals( tCounter++) = mMeshCluster->get_void_cells_in_cluster()( Ik )->get_index();
                }
            	return tVals;
		    }

//-----------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------

        };
    } /* namespace VIS */
} /* namespace moris */

#endif /* SRC_FEM_CL_VIS_VISUALIZATION_OBJECT_HPP_ */
