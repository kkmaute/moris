/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Factory.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_FACTORY_HPP_
#define SRC_HMR_CL_HMR_FACTORY_HPP_

#include "cl_HMR_Background_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src

namespace moris::hmr
{
    /**
     * \brief factory class that generates pointers to templated meshes
     */
    class Factory
    {
    private:

        // Stored parameters for building classes
        const Parameters* mParameters;

    public:

        /**
         * Constructor which takes in the HMR parameters needed to construct classes
         *
         * @param aParameters HMR parameters to use
         */
        explicit Factory( const Parameters* aParameters );

        /**
         * creates a background mesh depending on the number of dimensions set
         *
         * @return Background_Mesh_Base*   pointer to new background mesh
         */
        Background_Mesh_Base* create_background_mesh();

        /**
         * creates a Lagrange mesh depending on the number of dimensions set
         *
         * @param[in] aBackgroundMesh       pointer to background mesh
         * @param[in] aPolynomialDegree     degree of Lagrange mesh
         *
         * @return Mesh* pointer to new Lagrange mesh
         */
        Lagrange_Mesh_Base*
        create_lagrange_mesh(
                Background_Mesh_Base*       aBackgroundMesh,
                Cell< BSpline_Mesh_Base * > aBSplineMeshes,
                uint                        aActivationPattern,
                luint                       aPolynomialDegree );

        /**
         * creates a Lagrange mesh depending on the number of dimensions set
         *
         * @param[in] aBackgroundMesh       pointer to background mesh
         * @param[in] aPolynomialDegree     degree of Lagrange mesh
         *
         * @return Mesh* pointer to new Lagrange mesh
         */
        BSpline_Mesh_Base*
        create_bspline_mesh(
                Background_Mesh_Base* aBackgroundMesh,
                uint                  aActivationPattern,
                luint                 aPolynomialDegree );

    };
}
#endif /* SRC_HMR_CL_HMR_FACTORY_HPP_ */
