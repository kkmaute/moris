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
        MSI::Equation_Object * Element_Factory::create_interpolation_element(       Element_Type                        aElementType,
                                                                              const moris::Cell< const mtk::Cell * >  & aInterpolationCell,
                                                                                    moris::Cell< Node_Base* >         & aNodes,
                                                                                    Set                               * aSet )
        {
            return new fem::Interpolation_Element( aElementType, aInterpolationCell, aNodes, aSet );
        }

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

//------------------------------------------------------------------------------
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
                    MORIS_ERROR( false, "Element_Factory::create_element - Not a double sideset" );
                    break;
            }
            return tElement;
        }

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
