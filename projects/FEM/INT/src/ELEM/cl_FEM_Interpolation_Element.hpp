/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Interpolation_Element.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_ELEMENT_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_ELEMENT_HPP_

#include "assert.h"
#include <cmath>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"                     //MTK/src
#include "cl_MSI_Equation_Object.hpp"          //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                    //FEM/INT/src
#include "cl_FEM_Node.hpp"                     //FEM/INT/src
#include "cl_FEM_IWG.hpp"                      //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"       //FEM/INT/src
#include "cl_MTK_Integrator.hpp"               //MTK/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src

namespace moris::fem
{
    class Set;
    class Cluster;
    class Field_Interpolation_Manager;

    //------------------------------------------------------------------------------
    /**
     * \brief element class that communicates with the mesh interface
     */
    class Interpolation_Element : public MSI::Equation_Object
    {

      protected:
        // pointer to the mesh cluster
        // NOTE: first cluster in list is the one that gets built in XTK
        // NOTE: the subsequent clusters usually correspond to clusters created from each of the VIS meshes
        Vector< std::shared_ptr< fem::Cluster > > mFemCluster;

        // pointer to the leader and follower mesh interpolation cell
        const mtk::Cell* mLeaderInterpolationCell;
        const mtk::Cell* mFollowerInterpolationCell;

        // pointer to the fem set
        Set* mSet;

        // element type
        Element_Type mElementType;

        friend class Element_Bulk;
        friend class Element_Sideset;
        friend class Element_Double_Sideset;
        friend class Element_Time_Sideset;
        friend class Element_Time_Boundary;
        friend class Element;
        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        Interpolation_Element() {};

        /**
         * constructor
         * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
         * @param[ in ] aMeshCluster cluster pointer from mtk mesh
         * @param[ in ] aNodes       cell of node pointers
         * @param[ in ] aSet         a fem set
         */
        Interpolation_Element(
                const Element_Type                aElementType,
                const Vector< const mtk::Cell* >& aInterpolationCell,
                const Vector< Node_Base* >&       aNodes,
                Set*                              aSet );

        //------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Interpolation_Element() override {};

        //------------------------------------------------------------------------------
        /**
         * get ip cell
         */
        const mtk::Cell*
        get_ip_cell( mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderInterpolationCell;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerInterpolationCell;
                }
                default:
                    MORIS_ERROR( false, "aIsLeader can only be leader or follower" );
                    return mLeaderInterpolationCell;
            }
        }

        //------------------------------------------------------------------------------
        /**
         * set cluster
         * @param[ in ] aCluster     pointer to a fem cluster
         * @param[ in ] aFemMeshIndex   mesh index
         */
        void set_cluster(
                std::shared_ptr< fem::Cluster > aCluster,
                const uint                      aFemMeshIndex );

        //------------------------------------------------------------------------------
        /**
         * get cluster
         *
         * @param[ in ] aIndex   mesh index
         *
         * @ return   const reference to shared pointer of cluster
         */
        const std::shared_ptr< fem::Cluster >& get_cluster( const uint aIndex );

        //------------------------------------------------------------------------------

        size_t get_num_clusters() const
        {
            return mFemCluster.size();
        };

        //------------------------------------------------------------------------------
        /**
         * fill mat pdv assembly vector
         */
        void fill_mat_pdv_assembly_vector();

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian on cluster
         */
        void compute_jacobian() override;

        //------------------------------------------------------------------------------
        /**
         * compute the residual on cluster
         */
        void compute_residual() override;

        //------------------------------------------------------------------------------
        /**
         * compute the jacobian and the residual on cluster
         */
        void compute_jacobian_and_residual() override;

        //------------------------------------------------------------------------------
        /**
         * compute dRdp
         */
        void compute_dRdp() override;

        //------------------------------------------------------------------------------
        /**
         * compute dQIdp
         */
        void compute_dQIdp_explicit() override;

        //------------------------------------------------------------------------------
        /**
         * compute dRdp and dQIdp explicit
         */
        void compute_dQIdp_explicit_implicit() override;

        //------------------------------------------------------------------------------
        /**
         * compute dQIdp
         */
        void compute_dQIdp_implicit() override;

        //------------------------------------------------------------------------------
        /**
         * compute dQIdu
         */
        void compute_dQIdu() override;

        //------------------------------------------------------------------------------
        /**
         * compute the quantities of interest on cluster
         */
        void compute_QI() override;

        //------------------------------------------------------------------------------

        /**
         * compute the quantity of interest on cluster
         * @param[ in ] aFemMeshIndex  index for the mesh used
         * @param[ in ] aFieldType  an enum for computation/return type
         *                          GLOBAL, NODAL, ELEMENTAL_INT, ELEMENTAL_AVG
         */
        void
        compute_quantity_of_interest(
                const uint           aFemMeshIndex,
                enum vis::Field_Type aFieldType ) override;

        //------------------------------------------------------------------------------
        /**
         * Reorder adjoint solution vector to fit ordering used within elements
         *
         * @return reordered adjoint
         */
        Matrix< DDRMat > reorder_adjoint_pdofs();

        //------------------------------------------------------------------------------

        /**
         * @brief compute elemental coefficients and set FIs for evaluating IQIs
         */
        void setup_IQI_computation();

        //------------------------------------------------------------------------------

        void
        compute_nodal_QIs_standard(
                const uint           aMeshIndex,
                mtk::Leader_Follower aLeaderOrFollowerSide = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------

        void
        compute_nodal_QIs_double_sided( const uint aMeshIndex );

        //------------------------------------------------------------------------------

        void
        compute_nodal_QIs_nonconformal_side( uint const tMapElement );

        //------------------------------------------------------------------------------
        /**
         * FIXME get rid of mNodalWeakBCs
         * gets the nodal values for the weak BCs
         */
        Matrix< DDRMat >&
        get_weak_bcs() override
        {
            return mNodalWeakBCs;
        }

        //------------------------------------------------------------------------------

        void populate_fields(
                Vector< std::shared_ptr< fem::Field > >& aFields,
                Vector< std::string > const &            tFieldIQINames ) override;

        //------------------------------------------------------------------------------
        /**
         * set the field interpolators coefficients
         */
        void set_field_interpolators_coefficients();

        //------------------------------------------------------------------------------
        void
        create_prop_adv_assembly_data(
                Vector< Vector< Matrix< DDRMat > > >&                                         aAdvPropWeights,
                const Vector< sint >&                                                         aAdvIds,
                const Vector< Vector< std::shared_ptr< gen::Design_Extraction_Operator > > >& aExtractionOperators );

        //------------------------------------------------------------------------------
    };    // class FEM::Interpolation_Element

    //------------------------------------------------------------------------------
}    // namespace moris::fem

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_ELEMENT_HPP_ */
