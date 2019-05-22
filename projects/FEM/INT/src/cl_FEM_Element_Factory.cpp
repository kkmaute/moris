#include "assert.hpp"
#include "cl_FEM_Element_Factory.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Bulk.hpp"         //FEM/INT/src
#include "cl_FEM_Element_Sideset.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Time_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Cluster.hpp" //FEM/INT/src

#include "cl_MSI_Equation_Object.hpp" //FEM/MSI/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        Element_Factory::Element_Factory(){}

//------------------------------------------------------------------------------

        Element_Factory::~Element_Factory(){}

//------------------------------------------------------------------------------

//        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type         aElementType,
//                                                                mtk::Cell    const * aCell,
//                                                                Cell< Node_Base* > & aNodes,
//                                                                Set                * aElementBlock)
//        {
//            MSI::Equation_Object * tElement = nullptr;
//            tElement = new Cluster( aElementType, aCell, aNodes, aElementBlock );
//            return tElement;
//        }

        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type                aElementType,
                                                                mtk::Cell_Cluster const   * aCellCluster,
                                                                moris::Cell< Node_Base* > & aNodes,
                                                                Set                       * aSet )
        {
            MSI::Equation_Object * tCluster = nullptr;
            tCluster = new Cluster( aElementType, aCellCluster, aNodes, aSet );
            return tCluster;
        }

        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type                aElementType,
                                                                mtk::Side_Cluster const   * aSideCluster,
                                                                moris::Cell< Node_Base* > & aNodes,
                                                                Set                       * aSet )
        {
            MSI::Equation_Object * tCluster = nullptr;
            tCluster = new Cluster( aElementType, aSideCluster, aNodes, aSet );
            return tCluster;
        }

//------------------------------------------------------------------------------

        fem::Element * Element_Factory::create_element( Element_Type         aElementType,
                                                        mtk::Cell    const * aCell,
                                                        Set                * aElementBlock,
                                                        Cluster            * aCluster )
        {
            fem::Element * tElement = nullptr;

            switch( aElementType )
            {
                case ( fem::Element_Type::BULK ):
                    tElement = new Element_Bulk( aCell, aElementBlock, aCluster );
                    break;

                case ( fem::Element_Type::SIDESET ):
                    tElement = new Element_Sideset( aCell, aElementBlock, aCluster );
                    break;

                case ( fem::Element_Type::TIME_SIDESET ):
                    tElement = new Element_Time_Sideset( aCell, aElementBlock, aCluster );
                    break;

                default:
                    MORIS_ERROR( false, "Element_Factory::create_element - No element type specified" );
                    break;
            }
            return tElement;
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
