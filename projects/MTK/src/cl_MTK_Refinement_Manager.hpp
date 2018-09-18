/*
 * cl_MTK_Refinement_Manager.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: messe
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_REFINEMENT_MANAGER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_REFINEMENT_MANAGER_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src

namespace moris
{
    namespace mtk
    {
        class Field;
//------------------------------------------------------------------------------

        class Refinement_Manager
        {
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Refinement_Manager(){};

//------------------------------------------------------------------------------

            /**
             * trivial destructor
             */
            ~Refinement_Manager(){};

//------------------------------------------------------------------------------

            /**
             * returns indices of elements that are intersected
             *
             * @param[ out ] aVolumeCellIndices   indices of cells which are inside
             * @param[ out ] aVolumeCellIndices   indices of intersected cells
             * @param[ in  ] aScalarField      field that contains values
             * @param[ in  ] aLowerBound       lower bound of intersection ( default: 0.0 )
             * @param[ in  ] aUpperBound       upper bound of intersection ( default: 0.0 )
             */
            void
            find_volume_and_surface_cells(
                    Mat< moris_index > & aVolumeCellIndices,
                    Mat< moris_index > & aSurfaceCellIndices,
                    const Field        * aScalarField,
                    const real           aLowerBound=0.0,
                    const real           aUpperBound=0.0 );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_REFINEMENT_MANAGER_HPP_ */
