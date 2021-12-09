/*
 * cl_FEM_Interpolation_Element.hpp
 *
 *  Created on: Apr 20, 2019
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_ELEMENT_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_ELEMENT_HPP_

#include "assert.h"
#include <cmath>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"                  //MTK/src
#include "cl_MSI_Equation_Object.hpp"       //FEM/MSI/src
#include "cl_FEM_Enums.hpp"                 //FEM/INT/src
#include "cl_FEM_Node.hpp"                  //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_MTK_Integrator.hpp"            //MTK/src
#include "cl_FEM_Element_Factory.hpp"       //FEM/INT/src

namespace moris
{
    namespace fem
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
                moris::Cell< std::shared_ptr< fem::Cluster > > mFemCluster;

                // pointer to the master and slave mesh interpolation cell
                const mtk::Cell * mMasterInterpolationCell;
                const mtk::Cell * mSlaveInterpolationCell;

                // pointer to the fem set
                Set * mSet;

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
                Interpolation_Element(){};

                /**
                 * constructor
                 * @param[ in ] aElementType enum for element type (BULK, SIDESET, ...)
                 * @param[ in ] aMeshCluster cluster pointer from mtk mesh
                 * @param[ in ] aNodes       cell of node pointers
                 * @param[ in ] aSet         a fem set
                 */
                Interpolation_Element(
                        const Element_Type                       aElementType,
                        const moris::Cell< const mtk::Cell * > & aInterpolationCell,
                        const moris::Cell< Node_Base* >              & aNodes,
                        Set                                    * aSet );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~Interpolation_Element(){};

                //------------------------------------------------------------------------------
                /**
                 * get ip cell
                 */
                const mtk::Cell * get_ip_cell( mtk::Master_Slave aIsMaster )
                {
                    switch( aIsMaster )
                    {
                        case mtk::Master_Slave::MASTER :
                        {
                            return mMasterInterpolationCell;
                        }
                        case mtk::Master_Slave::SLAVE :
                        {
                            return mSlaveInterpolationCell;
                        }
                        default:
                            MORIS_ERROR( false, "aIsMaster can only be master or slave" );
                            return mMasterInterpolationCell;
                    }
                }

                //------------------------------------------------------------------------------
                /**
                 * set cluster
                 * @param[ in ] aCluster     pointer to a fem cluster
                 * @param[ in ] aMeshIndex   mesh index
                 */
                void set_cluster(
                        std::shared_ptr< fem::Cluster > aCluster,
                        const uint                      aMeshIndex );

                //------------------------------------------------------------------------------
                /**
                 * get cluster
                 *
                 * @param[ in ] aIndex   mesh index
                 *
                 * @ return   const reference to shared pointer of cluster
                 */
                const std::shared_ptr< fem::Cluster > & get_cluster( const uint aIndex );

                //------------------------------------------------------------------------------
                /**
                 * fill mat pdv assembly vector
                 */
                void fill_mat_pdv_assembly_vector();

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian on cluster
                 */
                void compute_jacobian();

                //------------------------------------------------------------------------------
                /**
                 * compute the residual on cluster
                 */
                void compute_residual();

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian and the residual on cluster
                 */
                void compute_jacobian_and_residual();

                //------------------------------------------------------------------------------
                /**
                 * compute dRdp
                 */
                void compute_dRdp();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp
                 */
                void compute_dQIdp_explicit();

                //------------------------------------------------------------------------------
                /**
                 * compute dRdp and dQIdp explicit
                 */
                void compute_dQIdp_explicit_implicit();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp
                 */
                void compute_dQIdp_implicit();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdu
                 */
                void compute_dQIdu();

                //------------------------------------------------------------------------------
                /**
                 * compute the quantities of interest on cluster
                 */
                void compute_QI();

                //------------------------------------------------------------------------------
                /**
                 * compute the quantity of interest on cluster
                 * @param[ in ] aMeshIndex  index for vis mesh used
                 * @param[ in ] aFieldType  an enum for computation/return type
                 *                          GLOBAL, NODAL, ELEMENTAL
                 */
                void compute_quantity_of_interest(
                        const uint                         aClusterIndex,
                        enum  vis::Field_Type              aFieldType );

                //------------------------------------------------------------------------------
                /**
                 * FIXME get rid of mNodalWeakBCs
                 * gets the nodal values for the weak BCs
                 */
                Matrix< DDRMat > & get_weak_bcs()
                {
                    return mNodalWeakBCs;
                }

                //------------------------------------------------------------------------------

                void populate_fields(
                        moris::Cell< std::shared_ptr< fem::Field > > & aFields,
                        moris::Cell< std::string > const             & tFieldIQINames );

                //------------------------------------------------------------------------------
            protected:

                //------------------------------------------------------------------------------
                /**
                 * set the field interpolators coefficients
                 */
                void set_field_interpolators_coefficients();

                //------------------------------------------------------------------------------
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTERPOLATION_ELEMENT_HPP_ */
