/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Equation_Object_Proxy.hpp
 *
 */

#pragma once

#include <memory>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_MTK_Enums.hpp"     //FEM/INT/src
#include "cl_MTK_Vertex.hpp"    //MTK/src

#include "fn_trans.hpp"
#include "op_times.hpp"

#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Equation_Object.hpp"
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
        class Equation_Object_Proxy : public Equation_Object
        {
            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------
            /**
             * trivial constructor
             */
            Equation_Object_Proxy() {};

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aElementBlock equation set pointer
             */
            Equation_Object_Proxy( Equation_Set* aEquationSet )
                    : Equation_Object( aEquationSet ) {};

            //------------------------------------------------------------------------------
            /**
             * constructor
             * @param[ in ] aNodeObjs leader/follower list of fem nodes
             */
            Equation_Object_Proxy( const Vector< Vector< fem::Node_Base* > >& aNodeObjs )
                    : Equation_Object( aNodeObjs )
            {
            }

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            virtual ~Equation_Object_Proxy() {};

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
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_jacobian - not implemented in msi." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * compute residual on equation object
             */
            virtual void
            compute_residual()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_residual - not implemented in msi." );
            }

            //-------------------------------------------------------------------------------------------------
            /**
             * compute jacobian and residual on equation object
             */
            virtual void
            compute_jacobian_and_residual()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_jacobian_and_residual - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dRdp on equation object
             */
            virtual void
            compute_dRdp()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_dRdp - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp explicit on equation object
             */
            virtual void
            compute_dQIdp_explicit()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_dQIdp_explicit - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp implicit on equation object
             */
            virtual void
            compute_dQIdp_implicit()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_dQIdp - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdp explicit and implicit on equation object
             */
            virtual void
            compute_dQIdp_explicit_implicit()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_dQIdp_explicit_implicit - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute dQIdu on equation object
             */
            virtual void
            compute_dQIdu()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_dQIdu - not implemented in msi." );
            }

            //------------------------------------------------------------------------------
            /**
             * compute QI on equation object
             */
            virtual void
            compute_QI()
            {
                MORIS_ERROR( false, "Equation_Object_Proxy::compute_QI - not implemented in msi." );
            };

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
                MORIS_ASSERT( false, "Equation_Object_Proxy::compute_quantity_of_interest() - not implemented for base class." );
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