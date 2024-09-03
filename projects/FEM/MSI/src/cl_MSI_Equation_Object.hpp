/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Object.hpp
 *
 */

#ifndef SRC_FEM_CL_EQUATION_OBJECT_HPP_
#define SRC_FEM_CL_EQUATION_OBJECT_HPP_

#include <memory>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Enums.hpp"     //FEM/INT/src
#include "cl_MTK_Vertex.hpp"    //MTK/src

#include "fn_trans.hpp"
#include "op_times.hpp"

#include "cl_MSI_Pdof_Host.hpp"
namespace moris
{
    class Dist_Vector;
    namespace mtk
    {
        class Set;
        class Cluster;
        class Field;
    }    // namespace mtk
    namespace fem
    {
        class Node_Base;
        class Element;
        class Field;
    }    // namespace fem
    namespace fem
    {
        class Cluster;
    }
    namespace vis
    {
        enum class Output_Type;
        enum class Field_Type;
    }    // namespace vis
    namespace MSI
    {
        struct Pdof;

        class Pdof_Host;
        class Pdv;
        class Pdv_Host;
        class Equation_Set;
        class Dof_Manager;
        class Equation_Object
        {
            //-------------------------------------------------------------------------------------------------

          protected:
            //-------------------------------------------------------------------------------------------------
            Vector< Vector< fem::Node_Base* > > mNodeObj;
            Vector< Vector< Pdof_Host* > >      mMyPdofHosts;    // Pointer to the pdof hosts of this equation object

            Vector< Pdof* >                     mFreePdofs;       // List of the pdof pointers of this equation obj
            Vector< Vector< Vector< Pdof* > > > mFreePdofList;    // FIXME list of free pdofs ordered after their dof type . mFreePdofs or mFreePdofList should be deleted

            Matrix< DDSMat >                     mUniqueAdofList;    // Unique adof list for this equation object
            Vector< Vector< Matrix< DDSMat > > > mUniqueAdofTypeList;
            moris::map< uint, uint >             mUniqueAdofMap;    // Map to

            Vector< Vector< moris::map< uint, uint > > > mUniqueAdofMapList;    // Map to

            //! weak BCs of element FIXME
            Matrix< DDRMat > mNodalWeakBCs;

            uint mEqnObjInd;

            Equation_Set* mEquationSet = nullptr;

            uint mNumPdofSystems = 0;

            // bool
            bool mUniqueAdofTypeListFlag = false;
            bool mFreePdofListFlag       = false;

            friend class fem::Element;
            friend class fem::Cluster;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Equation_Object() {};

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aElementBlock equation set pointer
             */
            Equation_Object( Equation_Set* aEquationSet )
                    : mEquationSet( aEquationSet ) {};

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aNodeObjs leader/follower list of fem nodes
             */
            Equation_Object( const Vector< Vector< fem::Node_Base* > >& aNodeObjs )
                    : mNodeObj( aNodeObjs )
            {
            }

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            virtual ~Equation_Object() {};

            //------------------------------------------------------------------------------
            /**
             * set time for equation object
             * @param[ in ] aTime matrix with time values to set
             */
            void
            set_time( Matrix< DDRMat >& aTime );

            //------------------------------------------------------------------------------
            /**
             * get time for equation object (from Equation model)
             * @param[ in ] mTime matrix with time values
             */
            Matrix< DDRMat >& get_time();

            //------------------------------------------------------------------------------
            /**
             * get previous time for equation object (from Equation model)
             * @param[ in ] mPreviousTime matrix with previous time values
             */
            Matrix< DDRMat >& get_previous_time();

            //------------------------------------------------------------------------------
            /**
             * @brief return the number of nodes, elements and ghosts related to this equation object.
             * This function is only for unit test purposes.
             */
            // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
            uint
            get_num_pdof_hosts()
            {
                uint tNumPdofHosts = 0;
                for ( uint Ik = 0; Ik < mNodeObj.size(); Ik++ )
                {
                    tNumPdofHosts = tNumPdofHosts + mNodeObj( Ik ).size();
                }
                return tNumPdofHosts;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief Returns the maximal pdof host (node) index of this equation object
             */
            uint get_max_pdof_hosts_ind();

            //------------------------------------------------------------------------------
            /**
             * create the pdof hosts for this equation object, if not created earlier
             * put the pdof hosts into the local pdof host list
             * This function is tested by the test [Eqn_Obj_create_pdof_host]
             * @param[in] aNumUsedDofTypes   Number of globally used dof types
             * @param[in] aPdofTypeMap       Map which maps the dof type enum values to
             *                               a consecutive list of dof type indices
             * @param[in] aPdofHostList      List of pdof hosts
             *
             */
            void create_my_pdof_hosts(
                    const uint              aNumUsedDofTypes,
                    const Matrix< DDSMat >& aPdofTypeMap,
                    const Matrix< DDUMat >& aTimePerDofType,
                    Vector< Pdof_Host* >&   aPdofHostList );

            //------------------------------------------------------------------------------

            const Vector< Vector< Pdof_Host* > >&
            get_pdof_hosts() const
            {
                return mMyPdofHosts;
            }

            //------------------------------------------------------------------------------
            /**
             * @brief create a list of pdof pointers related to this equation object
             * This function is tested by the test [Eqn_Obj_create_my_pdof_list]
             * [Dof_Mgn_create_unique_dof_type_map_matrix]
             */
            void create_my_pdof_list();

            //------------------------------------------------------------------------------
            /**
             * @brief create a unique list of adofs Ids corresponding to this equation object
             * This function is tested by the test [Eqn_Obj_create_my_list_of_adof_ids]
             */
            void create_my_list_of_adof_ids();

            //------------------------------------------------------------------------------
            /**
             * @brief create a map relating the adof ids to the positions for this equation object
             *  This function is tested by the test [Eqn_Obj_create_adof_map]
             */
            void set_unique_adof_map();

            //------------------------------------------------------------------------------
            /**
             * @brief create a PADofMap witch can be used to for a calculation from pdofs to adofs
             * This function is tested by the test [Eqn_Obj_PADofMap]
             */
            void build_PADofMap( Matrix< DDRMat >& aPADofMap );

            //------------------------------------------------------------------------------

            void build_PADofMap_list( Vector< Vector< Matrix< DDRMat > > >& aPADofMap );

            void build_PADofMap_1( Matrix< DDRMat >& aPADofMap );

            //------------------------------------------------------------------------------
            /**
             * @brief compute function for the pdof values of this particular equation object
             */
            void compute_my_pdof_values();

            //------------------------------------------------------------------------------
            /**
             * @brief compute function for the previous pdof values of this particular equation object
             */
            void compute_previous_pdof_values();

            //------------------------------------------------------------------------------
            /**
             * @brief compute function for the eigen vector pdof values of this particular equation object
             */
            void compute_my_eigen_vector_values();

            //------------------------------------------------------------------------------
            /**
             * @brief compute function for the adjoint values
             */
            void compute_my_adjoint_values();

            //------------------------------------------------------------------------------
            /**
             * @brief compute function for the previous adjoint values
             */
            void compute_my_previous_adjoint_values();

            //------------------------------------------------------------------------------
            /**
             * get the pdof values of this particular equation object.
             * get_my_pdof_values() has to be called first to initialize.
             * @param[ in ] aPdofValues           All pdof values of this equation object
             * @param[ in ] aRequestedDofTypes    List of requested dof types
             * @param[ in ] aRequestedPdofValues  Reference to the matrix of requested pdof values
             * @param[ in ] aIsLeader             enum for leader or follower
             */
            void get_my_pdof_values(
                    const Vector< Matrix< DDRMat > >&     aPdofValues,
                    const Vector< enum Dof_Type >&        aRequestedDofTypes,
                    Vector< Vector< Matrix< DDRMat > > >& aRequestedPdofValues,
                    const mtk::Leader_Follower            aIsLeader = mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            /**
             * reshape the pdof values of this equation object
             * @param[ in ] aPdofValues           list of matrices with pdof values
             *                                    (one matrix per dof type)
             * @param[ in ] aReshapedPdofValues   matrix with pdof values
             *                                    (one column per dof type)
             */
            void reshape_pdof_values(
                    const Vector< Matrix< DDRMat > >& aPdofValues,
                    Matrix< DDRMat >&                 aReshapedPdofValues );

            //------------------------------------------------------------------------------
            /**
             * reshape the pdof values of this equation object
             * @param[ in ] aPdofValues           list of matrices with pdof values
             *                                    (one matrix per dof type)
             * @param[ in ] aReshapedPdofValues   vector with pdof values
             */
            void reshape_pdof_values_vector(
                    const Vector< Matrix< DDRMat > >& aPdofValues,
                    Matrix< DDRMat >&                 aReshapedPdofValues );

            //------------------------------------------------------------------------------
            /**
             * FIXME doc????
             */
            void set_vector_entry_number_of_pdof();

            //------------------------------------------------------------------------------
            /**
             * get jacobian for equation object
             * @param[ in ] aEqnObjMatrix matrix to fill with jacobian on equation object
             */
            void get_egn_obj_jacobian( Matrix< DDRMat >& aEqnObjMatrix );

            //------------------------------------------------------------------------------
            /**
             * get residual on equation object
             * @param[ in ] aEqnObjRHS list of matrices to fill with RHS on equation object
             */
            void get_equation_obj_residual( Vector< Matrix< DDRMat > >& aEqnObjRHS );

            //------------------------------------------------------------------------------
            /**
             * get additional residual for staggered case on equation object
             * @param[ in ] aEqnObjRHS list of matrices to fill with RHS on equation object
             */
            void get_staggered_equation_obj_residual( Vector< Matrix< DDRMat > >& aEqnObjRHS );

            //------------------------------------------------------------------------------
            /**
             * get off-diagonal residual on equation object
             * @param[ in ] aEqnObjRHS list of matrices to fill with off-diagonal RHS on equation object
             */
            void get_equation_obj_off_diagonal_residual( Vector< Matrix< DDRMat > >& aEqnObjRHS );

            //-------------------------------------------------------------------------------------------------
            /**
             * get jacobian and residual on equation object
             * @param[ in ] aEqnObjMatrix matrix to fill with jacobian on equation object
             * @param[ in ] aEqnObjRHS list of matrices to fill with RHS on equation object
             */
            void get_egn_obj_jacobian_and_residual(
                    Matrix< DDRMat >&           aEqnObjMatrix,
                    Vector< Matrix< DDRMat > >& aEqnObjRHS );

            //-------------------------------------------------------------------------------------------------
            /**
             * add staggered contribution to residual
             * @param[ in ] aElementResidual ???
             */
            void add_staggered_contribution_to_residual( Vector< Matrix< DDRMat > >& aElementResidual );

            //-------------------------------------------------------------------------------------------------
            void get_equation_obj_dof_ids( Matrix< DDSMat >& aEqnObjAdofId );

            //-------------------------------------------------------------------------------------------------
            /**
             * returns a moris::Mat with indices of vertices that are connected to this element
             */
            moris_index get_node_index( const moris_index aElementLocalNodeIndex ) const;

            //-------------------------------------------------------------------------------------------------

            virtual Matrix< DDSMat >
            get_adof_indices()
            {
                MORIS_ERROR( false, "this function does nothing" );
                return Matrix< DDSMat >( 0, 0 );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * compute jacobian on equation object
             */
            virtual void
            compute_jacobian()
            {
                MORIS_ERROR( false, "Equation_Object::compute_jacobian - not implemented in msi." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * compute residual on equation object
             */
            virtual void
            compute_residual()
            {
                MORIS_ERROR( false, "Equation_Object::compute_residual - not implemented in msi." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * compute jacobian and residual on equation object
             */
            virtual void
            compute_jacobian_and_residual()
            {
                MORIS_ERROR( false, "Equation_Object::compute_jacobian_and_residual - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dRdp on equation object
             */
            virtual void
            compute_dRdp()
            {
                MORIS_ERROR( false, "Equation_Object::compute_dRdp - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp explicit on equation object
             */
            virtual void
            compute_dQIdp_explicit()
            {
                MORIS_ERROR( false, "Equation_Object::compute_dQIdp_explicit - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp implicit on equation object
             */
            virtual void
            compute_dQIdp_implicit()
            {
                MORIS_ERROR( false, "Equation_Object::compute_dQIdp - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp explicit and implicit on equation object
             */
            virtual void
            compute_dQIdp_explicit_implicit()
            {
                MORIS_ERROR( false, "Equation_Object::compute_dQIdp_explicit_implicit - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdu on equation object
             */
            virtual void
            compute_dQIdu()
            {
                MORIS_ERROR( false, "Equation_Object::compute_dQIdu - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute QI on equation object
             */
            virtual void
            compute_QI()
            {
                MORIS_ERROR( false, "Equation_Object::compute_QI - not implemented in msi." );
            };

            //------------------------------------------------------------------------------
            /**
             * compute integration error
             */
            virtual moris::real
            compute_integration_error(
                    moris::real ( *aFunction )( const Matrix< DDRMat >& aPoint ) )
            {
                MORIS_ERROR( false, "Equation_Object::compute_integration_error - not implemented in msi." );
                return 0.0;
            }

            //------------------------------------------------------------------------------
            /**
             * compute element average of scalar field
             */
            virtual moris::real
            compute_element_average_of_scalar_field()
            {
                MORIS_ERROR( false, "Equation_Object::compute_element_average_of_scalar_field - not implemented in msi." );
                return 0.0;
            }

            //------------------------------------------------------------------------------
            /**
             * return Neumann boundary conditions, writable version
             */
            virtual Matrix< DDRMat >&
            get_weak_bcs()
            {
                return mNodalWeakBCs;
            }

            //------------------------------------------------------------------------------
            /**
             * return Neumann boundary conditions, const version
             */
            const Matrix< DDRMat >&
            get_weak_bcs() const
            {
                return mNodalWeakBCs;
            }

            //------------------------------------------------------------------------------
            /**
             * return how many nodes are connected to this element
             */
            uint
            get_num_nodes() const
            {
                return mNodeObj( 0 ).size();
            }

            //------------------------------------------------------------------------------
            /**
             * get nodal pdof value
             * @param[ in ] aVertexIndex index for nodal value to get
             * @param[ in ] aDofType     list of dof type to get
             */
            virtual moris::real
            get_element_nodal_pdof_value(
                    moris_index                    aVertexIndex,
                    const Vector< MSI::Dof_Type >& aDofType )
            {
                MORIS_ERROR( false, "Equation_Object::get_element_nodal_pdof_value - this function does nothing" );
                return 0.0;
            }

            //------------------------------------------------------------------------------
            /**
             * set visualization cluster
             * @param[ in ] aVisMeshCluster mesh cluster pointer to set
             */
            virtual void
            set_visualization_cluster( const mtk::Cluster* aVisMeshCluster )
            {
                MORIS_ASSERT( false, "Equation_Object::set_visualization_cluster() - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute quantity of interest
             * @param[ in ] aFemMeshIndex mesh index to specify on which mesh to compute QI
             * @param[ in ] aFieldType enum for computation type (GLOBAL,NODAL,ELEMENTAL_INT,ELEMENTAL_AVG)
             */
            virtual void
            compute_quantity_of_interest(
                    const uint           aFemMeshIndex,
                    enum vis::Field_Type aFieldType )
            {
                MORIS_ASSERT( false, "Equation_Object::compute_quantity_of_interest() - not implemented for base class." );
            }

            //------------------------------------------------------------------------------

            virtual void
            populate_fields(
                    Vector< std::shared_ptr< fem::Field > >& aFields,
                    Vector< std::string > const &            aFieldIQINames )
            {
                MORIS_ASSERT( false, "Equation_Set::create_fields - not implemented for base class." );
            }

            //------------------------------------------------------------------------------
        };
    }    // namespace MSI
}    // namespace moris

#endif /* SRC_FEM_CL_EQUATION_OBJECT_HPP_ */
