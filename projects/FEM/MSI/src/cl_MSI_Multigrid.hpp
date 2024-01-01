/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Multigrid.hpp
 *
 */

#ifndef SRC_FEM_CL_MSI_MULTIGRID_HPP_
#define SRC_FEM_CL_MSI_MULTIGRID_HPP_

#include "cl_Vector.hpp"
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
    namespace dla
    {
        class Geometric_Multigrid;
    }
    namespace MSI
    {
        class Model_Solver_Interface;
        class Multigrid
        {
        private:
            //! Number of multigrid levels
            moris::uint mMultigridLevels = -1;

            //! Number of dofs which are on this level or coarser.
            moris::Matrix< DDUMat > mNumDofsRemain;

            //! List of external indices for each level
            Vector< Matrix< DDUMat > > mListAdofExtIndMap;

            //! List of type/time identifiers for each level
            Vector< Matrix< DDSMat > > mListAdofTypeTimeIdentifier;

            //! Map which maps external indices to internal MSI indices. List 1 = Level; List 2 = type/time;
            Vector< Vector< Matrix< DDSMat > > > mMultigridMap;

            // Mesh pointer
            mtk::Mesh * mMesh;

            // Pointer to the model solver interface
            moris::MSI::Model_Solver_Interface * mModelSolverInterface;

            // Maximal number of used tof types/time
            moris::sint mMaxDofTypes = -1;

        public:
            Multigrid( moris::MSI::Model_Solver_Interface * aModelSolverInterface,
                       moris::mtk::Mesh                   * aMesh );

            ~Multigrid(){};

            /**
             * @brief Initializing the member variable lists for the fines level.
             *
             */
            void multigrid_initialize();

            /**
             * @brief Create the member variable lists for the coarser levels
             *
             */
            void create_multigrid_level_dof_ordering();

            /**
             * @brief Create reverse map which maps external indices to the internal indices for each type/time and level.
             *
             */
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

            const Vector< Matrix< DDUMat > > & get_lists_of_ext_index_multigrid( )
            {
                return mListAdofExtIndMap;
            };

            const Vector< Matrix< DDSMat > > & get_lists_of_multigrid_identifiers( )
            {
                return mListAdofTypeTimeIdentifier;
            };

            const Vector< Vector< Matrix< DDSMat > > > & get_multigrid_map( )
            {
                return mMultigridMap;
            };

            const Matrix< DDUMat > & get_number_remaining_dofs( )
            {
                return mNumDofsRemain;
            };

        };
    } /* namespace MSI */
} /* namespace moris */

#endif /* SRC_FEM_CL_MSI_MULTIGRID_HPP_ */

