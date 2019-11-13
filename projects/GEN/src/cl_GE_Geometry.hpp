/*
 * cl_GE_Geometry.hpp
 *
 *  Created on: Dec 28, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GEOMETRY_HPP_
#define PROJECTS_GEN_SRC_CL_GEOMETRY_HPP_

// fem includes
#include "cl_FEM_Enums.hpp"

// ge includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Enums.hpp"
#include "cl_GE_Geometry_Library.hpp"

// hmr includes
#include "cl_HMR_Field.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"

// mtk includes
#include "cl_MTK_Mesh_Manager.hpp"

// other includes
#include <cmath>
#include <functional>
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
        virtual uint get_number_of_sub_types()
        {
            MORIS_ASSERT(false, "ge::Geometry::get_number_of_sub_types(): not implemented ");
            return 0;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief check that the analytic function has been set before attempting to use
         */
        virtual bool
        check_if_function_is_set( const moris_index aSubIndex )
        {
            MORIS_ASSERT(false, "ge::Geometry::check_if_function_is_set(): not implemented ");
            return false;
        };
        //------------------------------------------------------------------------------
        /*
         * @brief check that the analytic function for sensitivity has been set before attempting to use
         */
        virtual bool
        check_if_sensitivity_function_is_set( const moris_index aSubIndex )
        {
            MORIS_ASSERT(false,"ge::Geometry::check_if_sensitivity_function_is_set(): not implemented ");
            return false;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief sets the constants necessary for the specific geometry rep
         */
        virtual void
        set_constants( moris::Cell< real > aMyConstants, moris_index aSubIndex )
        {
            MORIS_ASSERT(false, "ge::Geometry::set_my_constants(): not implemented for this type of geometry representation");
        }
        //------------------------------------------------------------------------------
        virtual moris_index
        set_analytical_function_and_dphi_dp( std::function< moris::real ( const Matrix< DDRMat >    & aPoint,
                                                                          const moris::Cell< real >   aConstant ) > aFunc01,
                                             std::function< Matrix< DDRMat > ( const Matrix< DDRMat >    & aPoint,
                                                                               const moris::Cell< real >   aConstant ) > aFunc02,
                                             moris::Cell< real > aConstants )
        {
            MORIS_ASSERT(false,"ge::Geometry::set_analytical_function_and_dphi_dp(): please specify your own analytic function and derivative(s) ");
            return 0;
        };
        //------------------------------------------------------------------------------
		virtual moris_index
		set_analytical_function( real ( *mFuncAnalytic )( const Matrix< DDRMat > & aPoint, moris::Cell< real> aConstant ), moris::Cell< real > aConstants = {{ 0 }} )
		{
		    MORIS_ASSERT(false,"ge::Geometry::set_analytical_function(): please specify your own analytic function");
		    return 0;
		};
		//------------------------------------------------------------------------------
		virtual void
		set_analytical_function( real ( *funcPointer )( const Matrix< DDRMat >        & aCoordinates,
                                                              Cell<Cell<moris::real>> & aCenter,
                                                              Cell<moris::real>       & aRadius,
                                                              Cell<moris::real>       & aLength,
                                                              Cell<Cell<moris::real>> & aAxis ) )
        {
		    MORIS_ASSERT(false,"ge::Geometry::set_analytical_function(): please choose a valid function");
        };
        //------------------------------------------------------------------------------
        virtual moris_index
        set_analytical_function( AnalyticType aGeomType, moris::Cell< real > aConstants )
        {
            MORIS_ASSERT(false,"ge::Geometry::set_analytical_function(): please choose a valid function");
            return 0;
        };

        //------------------------------------------------------------------------------
        virtual void
        set_analytical_function_dphi_dp( Matrix< DDRMat > ( *mFuncAnalyticDphiDx )( const Matrix< DDRMat > & aPoint, Cell< real > aConst ), moris_index aFuncIndex = 0 )
        {
            MORIS_ASSERT(false,"ge::Geometry::set_analytical_function_dphi_dx(): please specify your own analytic function dphi/dx");
        };

        //------------------------------------------------------------------------------
        virtual void
        set_analytical_function_dphi_dp( AnalyticType aGeomType )
        {
            MORIS_ASSERT(false,"ge::Geometry::set_analytical_function_dphi_dx(): please choose a valid dphi/dx function");
        };

        //------------------------------------------------------------------------------
        virtual real
        get_field_val_at_coordinate( const Matrix< DDRMat > & aPoint,
                                     const moris_index aSubIndex = 0 )
        {
            MORIS_ASSERT(false,"ge::Geometry::get_field_val_at_coordinate(): function not implemented");
            return 0.0;
        };

        //------------------------------------------------------------------------------

        virtual Matrix< DDRMat >
        get_sensitivity_dphi_dp_at_coordinate( const Matrix< DDRMat >  & aPoint,
                                               const moris_index aSubIndex = 0 )
        {
            MORIS_ASSERT(false,"ge::Geometry::get_sensitivity_dphi_dp_at_coordinate(): function not implemented");
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
            MORIS_ASSERT(false,"ge::Geometry::get_node_adv_indices(): function not implemented");
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
            MORIS_ASSERT(false,"ge::Geometry::access_field_value_with_entity_index(): function not implemented");
            return 0.0;
        }

        //------------------------------------------------------------------------------
        virtual moris::Matrix< moris::DDRMat >
        evaluate_sensitivity_dphi_dp(moris::Matrix< moris::DDRMat > const & aLocalCoordinate, moris::size_t aEntityIndex, enum EntityRank aEntityRank)
        {
            MORIS_ASSERT(false,"ge::Geometry::evaluate_sensitivity_dx_dp(): function not implemented");

            moris::Matrix< moris::DDRMat > tSensitivityDxDp(1,1,0);
            return tSensitivityDxDp;
        }

        //------------------------------------------------------------------------------
        virtual std::string const &
        get_active_level_set_field_name() const
        {
            MORIS_ASSERT(false, "ge::Geometry::get_active_level_set_field_name(): function not implemented");

            return mDummyReturnString;
        }

        //------------------------------------------------------------------------------
        virtual void
        set_member_variables(moris::mtk::Mesh_Manager*   aMeshWithLevelSetFields,
                             Cell<std::string> const   & aFieldNames)
        {
            MORIS_ASSERT(false, "ge::Geometry::set_member_variables(): not implemented");
        }

        //------------------------------------------------------------------------------
        //fixme not sure what to pass in here, need generalized MTK field class or something similar
        //      to directly add the field to the mesh
        virtual void
        set_my_target_field( std::shared_ptr< hmr::Field > &aField )
        {
            MORIS_ASSERT(false, "ge::Geometry::set_my_target_field(): not implemented");
        }
        //------------------------------------------------------------------------------
        virtual void
        set_my_output_field( std::shared_ptr< hmr::Field > &aField )
        {
            MORIS_ASSERT(false, "ge::Geometry::set_my_target_field(): not implemented");
        }
        //------------------------------------------------------------------------------
        virtual std::shared_ptr< hmr::Field >
        get_my_target_field()
        {
            MORIS_ASSERT( false, "ge::Geometry::get_my_target_field(): not implemented" );
            return mDummyField;
        }
        //------------------------------------------------------------------------------
        virtual std::shared_ptr< hmr::Field >
        get_my_output_field()
        {
            MORIS_ASSERT( false, "ge::Geometry::get_my_output_field(): not implemented" );
            return mDummyField;
        }
        /*
         * *****************************************************************************
         * ************************* SDF GEOMETRY FUNCTIONS ****************************
         * *****************************************************************************
         */
        virtual moris_index
        add_hmr_field( std::shared_ptr< hmr::Field > &aField )
        {
            MORIS_ASSERT(false, "ge::Geometry::add_hmr_field(): not implemented");
            return 0;
        }
        //------------------------------------------------------------------------------
        virtual void
        initialize_sdf( const std::string & aObjectPath,
                        std::shared_ptr< mtk::Mesh > aMTKMesh,
                        const bool aVerboseFlag = true,
                        const uint aFieldIndex = 0 )
        {
            MORIS_ASSERT(false, "ge::Geometry::initalize_sdf(): not implemented");
        }
        //------------------------------------------------------------------------------
        virtual real
        get_sdf_vals( moris_index )
        {
            MORIS_ASSERT(false, "ge::Geometry::get_sdf_val_at_node_index(): not implemented");
            return 0.0;
        }

        /*
         * *****************************************************************************
         * ************************* FUNCTIONS FOR ALL TYPES ***************************
         * *****************************************************************************
         */
        //------------------------------------------------------------------------------
        /*
         * @brief set the mesh for the geometry representation
         */
        virtual void
        set_my_mesh(mtk::Mesh_Manager* aMyMesh)
        {
            MORIS_ASSERT(false, "ge::Geometry::set_mesh(): mesh not set");
        }
        //------------------------------------------------------------------------------
        /*
         * @brief sets the interpolation type and rule in both space and time, if these are not directly set, they are defaulted to linear Legrange in
         *        space and constant in time
         */
        virtual void
        set_my_interpolation_rules( fem::Interpolation_Type  aSpaceInterpType,
                                    mtk::Interpolation_Order aSpaceInterpOrder,
                                    fem::Interpolation_Type  aTimeInterpType  = fem::Interpolation_Type::CONSTANT,
                                    mtk::Interpolation_Order aTimeInterpOrder = mtk::Interpolation_Order::CONSTANT )
        {
            MORIS_ERROR( false, "ge::Geometry::set_my_interpolation_rules(): not implemented" );
        }
        //------------------------------------------------------------------------------
        /*
         * @brief add fields to the mesh associated with the geometry representation
         */
        virtual void
        add_field()
        {
            MORIS_ASSERT(false, "ge::Geometry::add_field(): not implemented");
        }
        //------------------------------------------------------------------------------
        /*
         * @brief function to report geometry representation type
         */
        virtual enum GeomType
        get_geom_type() const
        {
            MORIS_ASSERT(false, "ge::Geometry::get_geom_type(): not implemented");
            return GeomType::END_ENUM;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns a pointer to the geometry object's mesh
         */
        virtual mtk::Mesh_Manager*
        get_my_mesh()
        {
            MORIS_ASSERT(false, "ge::Geometry::get_my_mesh(): mesh has not been set");
            return mDummyMeshPointer;
        }
        //------------------------------------------------------------------------------
        /*
         * @brief returns the geometry object's specific constants
         */
        virtual moris::Cell< real >
        get_my_constants()
        {
            MORIS_ASSERT(false, "ge::Geometry::get_my_constants(): constants not set");
            return 0.0;
        }
        //------------------------------------------------------------------------------
        virtual fem::Interpolation_Type
        get_my_space_interpolation_type()
        {
            return fem::Interpolation_Type::LAGRANGE;
        }
        //------------------------------------------------------------------------------
        virtual mtk::Interpolation_Order
        get_my_space_interpolation_order()
        {
            return mtk::Interpolation_Order::LINEAR;
        }
        //------------------------------------------------------------------------------
        virtual fem::Interpolation_Type
        get_my_time_interpolation_type()
        {
            return fem::Interpolation_Type::LAGRANGE;
        }
        //------------------------------------------------------------------------------
        virtual mtk::Interpolation_Order
        get_my_time_interpolation_order()
        {
            return mtk::Interpolation_Order::LINEAR;
        }
	//------------------------------------------------------------------------------
    private:
        // dummy member variables
        std::string mDummyReturnString;
        Cell<std::string> mDummyReturnCell;
        mtk::Mesh_Manager* mDummyMeshPointer = nullptr;
        Matrix< DDRMat > mDummyMatrix;
        std::shared_ptr< hmr::Field > mDummyField = nullptr;
	//------------------------------------------------------------------------------
    protected:

	};


//------------------------------------------------------------------------------

} /* namespace gen */
} /* namespace moris */


#endif /* PROJECTS_GEN_SRC_CL_GEOMETRY_HPP_ */
