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

        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type         aElementType,
                                                                mtk::Cell    const * aCell,
                                                                Cell< Node_Base* > & aNodes,
                                                                Element_Block      * aElementBlock)
        {
            MSI::Equation_Object * tElement = nullptr;

//            switch( aElementType )
//            {
//                case ( fem::Element_Type::BULK ):
                    tElement = new Cluster( aElementType, aCell, aNodes, aElementBlock );
//                    break;

//                case ( fem::Element_Type::BULK ):
//                    tElement = new Element_Bulk( aCell, aNodes, aElementBlock );
//                    break;
//
//                case ( fem::Element_Type::SIDESET ):
//                    tElement = new Element_Sideset( aCell, aNodes, aElementBlock );
//                    break;
//
//                case ( fem::Element_Type::TIME_SIDESET ):
//                    tElement = new Element_Time_Sideset( aCell, aNodes, aElementBlock );
//                    break;

//                default:
//                    MORIS_ERROR( false, "Element_Factory::create_cluster - No element type specified" );
//                    break;
//            }
            return tElement;
        }

//------------------------------------------------------------------------------

        fem::Element * Element_Factory::create_element( Element_Type         aElementType,
                                                        mtk::Cell    const * aCell,
                                                        Element_Block      * aElementBlock,
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
