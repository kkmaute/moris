/*
 * cl_MSI_Multigrid.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_MSI_MULTIGRID_HPP_
#define SRC_FEM_CL_MSI_MULTIGRID_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#include "cl_Map.hpp"
#include "fn_sum.hpp"

namespace moris
{
    namespace mtk
    {
        class Mesh;
    }
    namespace MSI
    {
        class Model_Solver_Interface;
        class Multigrid
        {
        private:
            moris::sint mMultigridLevels = -1;

            mtk::Mesh * mMesh;

            moris::MSI::Model_Solver_Interface * mModelSolverInterface;

            moris::sint mMaxDofTypes = -1;

            moris::Cell< Matrix< DDUMat > > mListAdofExtIndMap;            // List of fine of coarse external index
            moris::Cell< Matrix< DDSMat > > mListAdofTypeTimeIdentifier;   // List of type time identifiers for coarse and fine mesh

            moris::Cell< moris::Cell< Matrix< DDSMat > > > mMultigridMap;  // Map which maps external indices to internal MSI indices. List 1 = Level; List 2 = type/time;


        public:
            Multigrid( moris::MSI::Model_Solver_Interface * aModelSolverInterface,
                       moris::mtk::Mesh                   * aMesh );

            ~Multigrid(){};

            void multigrid_initialize();

            void create_multigrid_level_dof_ordering();

            void create_multigrid_maps();

            /**
             * @brief Function to read internal numbering from multigrid maps based on external indices
             *
             * @param[in] aLevel                 Level for which information accounts. Fines mesh in the V-Cycle is 0
             * @param[in] aExtFineIndices        List with external indices.
             * @param[in] aTypeTimeIdentifier    Type and time identifier.
             *
             * @param[out] aInternalFineIndices  Internal indices for this type and time corresponding to inputs.
             *
             */
            void read_multigrid_maps( const moris::uint               aLevel,
                                      const moris::Matrix< DDSMat > & aExtFineIndices,
                                      const moris::sint               aTypeTimeIdentifier,
                                            moris::Matrix< DDSMat > & aInternalFineIndices);

            moris::Cell< Matrix< DDUMat > > get_lists_of_ext_index_multigrid( )
            {
                return mListAdofExtIndMap;
            };

        };
    } /* namespace MSI */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_MULTIGRID_HPP_ */
