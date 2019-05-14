/*
 * cl_GE_Geometry.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GEOMETRY_HPP_
#define PROJECTS_GEN_SRC_CL_GEOMETRY_HPP_

// ge includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Enums.hpp"
#include "cl_GE_Interface.hpp"
#include "cl_GE_Geometry_Library.hpp"
//------------------------------------------------------------------------------
// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"
//------------------------------------------------------------------------------
// mtk includes
#include "cl_MTK_Mesh_Manager.hpp"
//------------------------------------------------------------------------------
// other includes
#include <cmath>
//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
	class Geometry
	{

	public:
		Geometry( )
	    {
	    };

		~Geometry(){};
		/*
		 * *****************************************************************************
		 * ************************ ANALYTIC GEOMETRY FUNCTIONS ************************
		 * *****************************************************************************
		 */
        //------------------------------------------------------------------------------
		virtual
		void set_analytical_function( real ( *mFuncAnalytic )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant ) )
		{
		    MORIS_ASSERT(false,"set_analytical_function(): please specify your own analytic function");
		};
		//------------------------------------------------------------------------------
		virtual
        void set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat >        & aCoordinates,
                                                                   Cell<Cell<moris::real>> & aCenter,
                                                                   Cell<moris::real>       & aRadius,
                                                                   Cell<moris::real>       & aLength,
                                                                   Cell<Cell<moris::real>> & aAxis ) )
        {
		    MORIS_ASSERT(false,"set_analytical_function(): please choose a valid function");
        };
        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function( type aGeomType )
        {
            MORIS_ASSERT(false,"set_analytical_function(): please choose a valid function");
        };

        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function_dphi_dx( Matrix< DDRMat > ( *mFuncAnalyticDphiDx )( const Matrix< DDRMat > & aPoint, Cell< real > aConst ) )
        {
            MORIS_ASSERT(false,"set_analytical_function_dphi_dx(): please specify your own analytic function dphi/dx");
        };

        //------------------------------------------------------------------------------
        virtual
        void set_analytical_function_dphi_dx( type aGeomType )
        {
            MORIS_ASSERT(false,"set_analytical_function_dphi_dx(): please choose a valid dphi/dx function");
        };

        //------------------------------------------------------------------------------
        virtual
        real get_field_val_at_coordinate( const Matrix< DDRMat > & aPoint, moris::Cell< real > aConst = 0.0 )
        {
            MORIS_ASSERT(false,"get_field_val_at_coordinate(): function not implemented");
            return 0.0;
        };

        //------------------------------------------------------------------------------

        virtual
        Matrix< DDRMat >
        get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint,
                                                     moris::Cell< real > aConst)
        {
            MORIS_ASSERT(false,"get_sensitivity_dphi_dp_at_coordinate(): function not implemented");
            Matrix< DDRMat > tSensitivityDxDp(4, 3, 0.0);
            return tSensitivityDxDp;
        };

        //------------------------------------------------------------------------------
        /*
         * *****************************************************************************
         * ************************ DISCRETE GEOMETRY FUNCTIONS ************************
         * *****************************************************************************
         */
        //------------------------------------------------------------------------------
        virtual moris::Matrix< moris::IndexMat >
        get_node_adv_indices(moris::Matrix< moris::IndexMat > const & aNodeIndices)
        {
            MORIS_ASSERT(false,"get_node_adv_indices(): function not implemented");
            Matrix< IndexMat > tTemp(1,1);
            tTemp(0,0) = 0.0;
            return tTemp;
        }

        //------------------------------------------------------------------------------
        /*
         * @brief returns the field value of a specific entity
         *
         * @param[in] aEntityIndex - index of the entity
         * @param[in] aEntityRank  - entity type (e.g. NODE)
         *
         * @param[out] field value
         */
        virtual moris::real
        access_field_value_with_entity_index(moris::moris_index aEntityIndex,
                                             enum EntityRank    aEntityRank) const
        {
            MORIS_ASSERT(false,"access_field_value_with_entity_index(): function not implemented");
            return 0.0;
        }

        //------------------------------------------------------------------------------
        virtual moris::Matrix< moris::DDRMat >
        evaluate_sensitivity_dphi_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum EntityRank aEntityRank)
        {
            MORIS_ASSERT(false,"evaluate_sensitivity_dx_dp(): function not implemented");

            moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
            return tSensitivityDxDp;
        }

        //------------------------------------------------------------------------------
        virtual std::string const &
        get_active_level_set_field_name() const
        {
            MORIS_ASSERT(false, "get_active_level_set_field_name(): function not implemented");

            return mDummyReturnString;
        }

        //------------------------------------------------------------------------------
        virtual moris::size_t
        get_num_levelset() const
        {
            MORIS_ASSERT(false, "get_num_levelset(): function not implemented");

            size_t tReturn = 1.0;
            return tReturn;
        }

        //------------------------------------------------------------------------------
        virtual bool
        advance_to_next_level_set()
        {
            MORIS_ASSERT(false, "advance_to_next_level_set(): function not implemented");

            return false;
        }

        //------------------------------------------------------------------------------
        virtual std::string const &
        get_level_set_field_name(moris::size_t aLevelSetIndex) const
        {
            MORIS_ASSERT(false, "get_level_set_field_name(): not implemented");

            return mDummyReturnString;
        }

        //------------------------------------------------------------------------------
        virtual Cell<std::string> const &
        get_level_set_field_name() const
        {
            MORIS_ASSERT(false, "get_level_set_field_name(): not implemented");

            return mDummyReturnCell;
        }


        //------------------------------------------------------------------------------
        virtual moris::mtk::Mesh*
        get_level_set_mesh()
        {
            MORIS_ASSERT(false, "get_level_set_mesh(): not implemented");

            return mDummyMeshPointer;
        }

        //------------------------------------------------------------------------------
        virtual void
        set_member_variables(moris::mtk::Mesh*         aMeshWithLevelSetFields,
                             Cell<std::string> const & aFieldNames)
        {
            MORIS_ASSERT(false, "set_member_variables(): not implemented");
        }
        /*
         * *****************************************************************************
         * ************************* FUNCTIONS FOR BOTH TYPES **************************
         * *****************************************************************************
         */
        //------------------------------------------------------------------------------
        virtual bool
        is_analytic() const
        {
//            MORIS_ASSERT(false, "is_analytic(): child geometry type not specified");
            return true;
        }

        //------------------------------------------------------------------------------
        /*
         * @brief set the mesh and T-matrix for the geometry representation
         */
        void
        set_mesh_and_t_matrix(mtk::Mesh_Manager & aMyMesh,
                              Matrix< DDRMat >  & aMyTMatrix)
        {
            mMyMesh    = & aMyMesh;
            mMyTMatrix = aMyTMatrix;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns a pointer to the geometry object's mesh
         */
        mtk::Mesh_Manager*
        get_my_mesh()
        {
            return mMyMesh;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns a pointer to the geometry object's T-matrix
         */
        Matrix< DDRMat >
        get_my_t_matrix()
        {
            return mMyTMatrix;
        }


	//------------------------------------------------------------------------------
    private:
        mtk::Mesh_Manager* mMyMesh;
        Matrix< DDRMat >   mMyTMatrix;

        std::string mDummyReturnString;
        Cell<std::string> mDummyReturnCell;
        mtk::Mesh* mDummyMeshPointer = nullptr;
        Matrix< DDRMat > mDummyMatrix;
	//------------------------------------------------------------------------------
    protected:

	};


//------------------------------------------------------------------------------

} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_GEN_SRC_CL_GEOMETRY_HPP_ */
