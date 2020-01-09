/*
 * cl_FEM_Element_Factory.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_

#include "assert.h"

#include "typedefs.hpp"               //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Enums.hpp"           //FEM/INT/src
#include "cl_FEM_IWG.hpp"             //FEM/INT/src
#include "cl_FEM_Node.hpp"            //FEM/INT/src
#include "cl_MTK_Cell.hpp" //MTK
#include "cl_MTK_Cell_Cluster.hpp" //MTK
#include "cl_MTK_Side_Cluster.hpp" //MTK
#include "cl_MTK_Double_Side_Cluster.hpp" //MTK

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
//        /**
//         * create a cell cluster
//         * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
//         * @param[ in ] aMeshCluster pointer to corresponding mesh cluster
//         * @param[ in ] aNodes       list of fem node pointers for IP vertices
//         * @param[ in ] aSet         pointer to corresponding fem set
//         */
//        MSI::Equation_Object * create_cluster( Element_Type                aElementType,
//                                               mtk::Cluster const        * aMeshCluster,
//                                               moris::Cell< Node_Base* > & aNodes,
//                                               Set                       * aSet );
        /**
         * create a cell cluster
         * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
         * @param[ in ] aMeshCluster pointer to corresponding mesh cluster
         * @param[ in ] aNodes       list of fem node pointers for IP vertices
         * @param[ in ] aIGNodes     list of fem node pointers for IG vertices
         * @param[ in ] aSet         pointer to corresponding fem set
         */
        MSI::Equation_Object * create_cluster( Element_Type                aElementType,
                                               mtk::Cluster const        * aMeshCluster,
                                               moris::Cell< Node_Base* > & aNodes,
                                               moris::Cell< Node_Base* > & aIGnodes,
                                               Set                       * aSet );

//------------------------------------------------------------------------------
        /**
         * create element
         */
        fem::Element * create_element( Element_Type         aElementType,
                                       mtk::Cell    const * aCell,
                                       Set                * aSet,
                                       Cluster            * aCluster,
                                       moris::moris_index   aCellIndexInCluster );

        fem::Element * create_element( Element_Type         aElementType,
                                       mtk::Cell    const * aLeftCell,
                                       mtk::Cell    const * aRightCell,
                                       Set                * aSet,
                                       Cluster            * aCluster,
                                       moris::moris_index   aCellIndexInCluster );

    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_FACTORY_HPP_ */
