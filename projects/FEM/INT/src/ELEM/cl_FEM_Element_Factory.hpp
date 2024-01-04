/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Factory.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_

#include "assert.h"

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Enums.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG.hpp"                    //FEM/INT/src
#include "cl_FEM_Node.hpp"                   //FEM/INT/src
#include "cl_MTK_Cell.hpp"                   //MTK
#include "cl_MTK_Cell_Cluster.hpp"           //MTK
#include "cl_MTK_Side_Cluster.hpp"           //MTK
#include "cl_MTK_Double_Side_Cluster.hpp"    //MTK

namespace moris
{
    namespace MSI
    {
        class Equation_Object;
    }
    //------------------------------------------------------------------------------
    namespace fem
    {
        class Set;
        class Element;
        class Cluster;
        //------------------------------------------------------------------------------

        /**
         * \brief element factory
         */
        class Element_Factory
        {

            //------------------------------------------------------------------------------
          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            Element_Factory();

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~Element_Factory();

            //------------------------------------------------------------------------------
            /**
             * create an interpolation element
             * @param[ in ] aElementType       enum for element type (BULK, SIDESET, ...)
             * @param[ in ] aInterpolationCell pointer to corresponding IP mesh cell
             * @param[ in ] aNodes             list of fem node pointers for IP vertices
             * @param[ in ] aSet               pointer to corresponding fem set
             */
            MSI::Equation_Object* create_interpolation_element(
                    Element_Type                      aElementType,
                    const Vector< const mtk::Cell* >& aInterpolationCell,
                    Vector< Node_Base* >&             aNodes,
                    Set*                              aSet );

            //------------------------------------------------------------------------------
            /**
             * create an integration element
             * @param[ in ] aElementType        enum for element type (BULK, SIDESET, ...)
             * @param[ in ] aCell               pointer to corresponding IG mesh cell
             * @param[ in ] aSet                pointer to corresponding fem set
             * @param[ in ] aCluster            pointer to corresponding fem cluster
             * @param[ in ] aCellIndexInCluster a index for the IG cell within the cluster
             */
            fem::Element* create_single_sided_element(
                    Element_Type       aElementType,
                    mtk::Cell const *  aCell,
                    Set*               aSet,
                    Cluster*           aCluster,
                    moris::moris_index aCellIndexInCluster );

            //------------------------------------------------------------------------------
            /**
             * create an integration element
             * @param[ in ] aElementType        enum for element type (BULK, SIDESET, ...)
             * @param[ in ] aLeftCell           pointer to corresponding leader IG mesh cell
             * @param[ in ] aRightCell          pointer to corresponding follower IG mesh cell
             * @param[ in ] aSet                pointer to corresponding fem set
             * @param[ in ] aCluster            pointer to corresponding fem cluster
             * @param[ in ] aCellIndexInCluster a index for the IG cell within the cluster
             */
            fem::Element* create_double_sided_element(
                    Element_Type       aElementType,
                    mtk::Cell const *  aLeftCell,
                    mtk::Cell const *  aRightCell,
                    Set*               aSet,
                    Cluster*           aCluster,
                    moris::moris_index aCellIndexInCluster );

        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_ */
