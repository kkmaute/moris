/*
 * cl_FEM_Element_Block.hpp
 *
 *  Created on: Mar 10, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_

#include "assert.h"
#include "cl_MSI_Equation_Block.hpp"               //FEM/INT/src
#include "cl_FEM_Enums.hpp"               //FEM/INT/src
#include "cl_MTK_Enums.hpp"               //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"               //FEM/INT/src

#include "cl_Communication_Tools.hpp"               //FEM/INT/src


namespace moris
{
namespace mtk
{
   class Cell;
}
namespace MSI
{
   class Model_Solver_Interface;
}
    namespace fem
    {
    class IWG;
    class Field_Interpolator;
    class Geometry_Interpolator;

//------------------------------------------------------------------------------
    /**
     * \brief element block class that communicates with the mesh interface
     */
    class Element_Block : public MSI::Equation_Block
    {
    private:
        moris::Cell< mtk::Cell const* >     mMeshElementPointer;

        moris::Cell< Node_Base* >           mNodes;

        Geometry_Interpolator             * mGeometryInterpolator = nullptr;

        moris::Cell< Field_Interpolator* >  mFieldInterpolators;

        // cell of pointers to IWG objects
        moris::Cell< IWG* > mIWGs;

        enum fem::Element_Type mElementType;

        // map of the element active dof types
        moris::Matrix< DDSMat >                   mInterpDofTypeMap;
        moris::Cell< Cell< enum MSI::Dof_Type > > mInterpDofTypeList;
        uint                                      mNumOfInterp;

        enum fem::Integration_Order mIntegrationOrder;

        enum fem::Integration_Order mSideIntegrationOrder;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * constructor
         *
         * @param[ in ]     List of mtk::Cell pointer
         */
        Element_Block( moris::Cell< mtk::Cell const * > & aCell,
                       Element_Type                aElementType,
                       Cell< IWG* >              & aIWGs,
                       Cell< Node_Base* >        & aNodes);

        /**
         * trivial constructor
         */
        Element_Block( ){};

//------------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Element_Block();

//------------------------------------------------------------------------------

        void delete_pointers();

//------------------------------------------------------------------------------

        void finalize( MSI::Model_Solver_Interface * aModelSolverInterface );

//------------------------------------------------------------------------------

        void create_dof_type_lists();

//------------------------------------------------------------------------------

        void create_unique_dof_type_lists();

//------------------------------------------------------------------------------

        void create_unique_list_of_first_dof_type_of_group();

//------------------------------------------------------------------------------

        Cell< MSI::Equation_Object * > & get_equation_object_list()
        {
            return mElements;
        };

//------------------------------------------------------------------------------

        moris::Cell< Field_Interpolator* > & get_block_field_interpolator()
        {
            return mFieldInterpolators;
        }

//------------------------------------------------------------------------------

        Geometry_Interpolator * get_block_geometry_interpolator()
        {
            return mGeometryInterpolator;
        }

//------------------------------------------------------------------------------

        uint get_num_IWG()
        {
            return mIWGs.size();
        }

//------------------------------------------------------------------------------

        moris::Cell< IWG* > & get_IWGs()
        {
            return mIWGs;
        }

//------------------------------------------------------------------------------

        moris::Cell< enum MSI::Dof_Type > & get_unique_dof_type_list()
        {
            return mEqnObjDofTypeList;
        }

//------------------------------------------------------------------------------

        moris::Matrix< DDSMat > & get_interpolator_dof_type_map()
        {
            return mInterpDofTypeMap;
        }

//------------------------------------------------------------------------------

        moris::Cell< Cell< enum MSI::Dof_Type > > & get_interpolator_dof_type_list()
        {
            return mInterpDofTypeList;
        }

//------------------------------------------------------------------------------

        uint & get_num_interpolators()
        {
            return mNumOfInterp;
        }

//------------------------------------------------------------------------------

        enum fem::Integration_Order & get_integration_order()
        {
            return mIntegrationOrder;
        }

//------------------------------------------------------------------------------

        enum fem::Integration_Order & get_side_integration_order()
        {
            return mSideIntegrationOrder;
        }

//------------------------------------------------------------------------------

        /**
         * get the field interpolators for an IWG
         */
        moris::Cell< Field_Interpolator* > get_IWG_field_interpolators ( IWG*                               & aIWG,
                                                                         moris::Cell< Field_Interpolator* > & aFieldInterpolators );

//------------------------------------------------------------------------------

        /**
         * auto detect full integration scheme
         */
        //FIXME: works for Lagrange only
        mtk::Interpolation_Order get_auto_interpolation_order( const moris::uint aNumVertices,
                                                               const mtk::Geometry_Type aGeometryType );

//------------------------------------------------------------------------------

        fem::Interpolation_Type get_auto_time_interpolation_type( const moris::uint aNumVertices );

//------------------------------------------------------------------------------

        void create_field_interpolators( MSI::Model_Solver_Interface * aModelSolverInterface );

        /**
          * auto detect interpolation scheme
          */
        fem::Integration_Order get_auto_integration_order( const mtk::Geometry_Type aGeometryType );

    };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_BLOCK_HPP_ */
