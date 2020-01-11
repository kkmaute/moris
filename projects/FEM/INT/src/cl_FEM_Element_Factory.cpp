#include "assert.hpp"
#include "cl_FEM_Element_Factory.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Bulk.hpp"         //FEM/INT/src
#include "cl_FEM_Element_Sideset.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Double_Sideset.hpp"      //FEM/INT/src
#include "cl_FEM_Element_Time_Sideset.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Element.hpp" //FEM/INT/src

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
        MSI::Equation_Object * Element_Factory::create_interpolation_element(       Element_Type                aElementType,
                                                                              const moris::Cell< const mtk::Cell * >  & aInterpolationCell,
                                                                                    moris::Cell< Node_Base* > & aNodes,
                                                                                    Set                       * aSet )
        {
            return new fem::Interpolation_Element( aElementType, aInterpolationCell, aNodes, aSet );
        }

//        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type                aElementType,
//                                                                mtk::Cell_Cluster const   * aCellCluster,
//                                                                moris::Cell< Node_Base* > & aNodes,
//                                                                Set                       * aSet )
//        {
//            MSI::Equation_Object * tCluster = nullptr;
//            tCluster = new Cluster( aElementType, aCellCluster, aNodes, aSet );
//            return tCluster;
//        }
//
//        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type                aElementType,
//                                                                mtk::Side_Cluster const   * aSideCluster,
//                                                                moris::Cell< Node_Base* > & aNodes,
//                                                                Set                       * aSet )
//        {
//            MSI::Equation_Object * tCluster = nullptr;
//            tCluster = new Cluster( aElementType, aSideCluster, aNodes, aSet );
//            return tCluster;
//        }
//
//        MSI::Equation_Object * Element_Factory::create_cluster( Element_Type                     aElementType,
//                                                                mtk::Double_Side_Cluster const & aDoubleSideCluster,
//                                                                moris::Cell< Node_Base* >      & aNodes,
//                                                                Set                            * aSet )
//        {
//            MSI::Equation_Object * tCluster = nullptr;
//            tCluster = new Cluster( aElementType, aDoubleSideCluster, aNodes, aSet );
//            return tCluster;
//        }

//------------------------------------------------------------------------------

        fem::Element * Element_Factory::create_element( Element_Type         aElementType,
                                                        mtk::Cell    const * aCell,
                                                        Set                * aSet,
                                                        Cluster            * aCluster,
                                                        moris::moris_index   aCellIndexInCluster )
        {
            fem::Element * tElement = nullptr;

            switch( aElementType )
            {
                case ( fem::Element_Type::BULK ):
                    tElement = new Element_Bulk( aCell, aSet, aCluster, aCellIndexInCluster );
                    break;

                case ( fem::Element_Type::SIDESET ):
                    tElement = new Element_Sideset( aCell, aSet, aCluster, aCellIndexInCluster );
                    break;

                case ( fem::Element_Type::TIME_SIDESET ):
                    tElement = new Element_Time_Sideset( aCell, aSet, aCluster, aCellIndexInCluster );
                    break;

                default:
                    MORIS_ERROR( false, "Element_Factory::create_element - No element type specified" );
                    break;
            }
            return tElement;
        }

        fem::Element * Element_Factory::create_element( Element_Type         aElementType,
                                                        mtk::Cell    const * aLeftCell,
                                                        mtk::Cell    const * aRightCell,
                                                        Set                * aSet,
                                                        Cluster            * aCluster,
                                                        moris::moris_index   aCellIndexInCluster )
        {
            fem::Element * tElement = nullptr;

            switch( aElementType )
            {
                case ( fem::Element_Type::DOUBLE_SIDESET ):
                    tElement = new Element_Double_Sideset( aLeftCell, aRightCell, aSet, aCluster, aCellIndexInCluster );
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
