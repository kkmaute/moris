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
#include "cl_FEM_Integrator.hpp"            //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"       //FEM/INT/src
#include "cl_FEM_Set.hpp"                   //FEM/INT/src


namespace moris
{
    namespace fem
    {
    class Set;
    class Cluster;
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
        Interpolation_Element( const Element_Type                       aElementType,
                               const moris::Cell< const mtk::Cell * > & aInterpolationCell,
                                     moris::Cell< Node_Base* >        & aNodes,
                                     Set                              * aSet );

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        ~Interpolation_Element(){};

//------------------------------------------------------------------------------
        /**
         * set cluster
         * @param[ in ] aCluster pointer to a fem cluster
         * @param[ in ] aIndex   mesh index
         */
        void set_cluster( std::shared_ptr< fem::Cluster > aCluster,
                          const uint                      aMeshIndex )
        {
            // if mesh index is 0 (i.e., forward analysis mesh, IG mesh)
            if( aMeshIndex == 0 )
            {
                // fem cluster with index 0 should be set only once and shall not be changed
                MORIS_ASSERT( !( mFemCluster.size() >= 1 ),
                              "Interpolation_Element::set_cluster() - first fem cluster is already set");
            }

            // get max size for fem cluster list
            sint tSize = std::max( ( sint )mFemCluster.size(), ( sint )aMeshIndex + 1 );

            // resize fem cluster list
            mFemCluster.resize( tSize );

            // add the fem cluster to the list
            mFemCluster( aMeshIndex ) = aCluster;
        };

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
         * compute the quantity of interest on cluster
         * @param[ in ] aMeshIndex  index for vis mesh used
         * @param[ in ] aOutputType an enum for output type
         * @param[ in ] aFieldType  an enum for computation/return type
         *                          GLOBAL, NODAL, ELEMENTAL
         */
        void compute_quantity_of_interest( const uint             aClusterIndex,
                                           enum  vis::Output_Type aOutputType,
                                           enum  vis::Field_Type  aFieldType );

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
        /**
         * get the nodal pdof values on an element
         * @param[ in ] aVertexIndex a vertex index
         * @param[ in ] aDofType     a dof type
         */
        real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                                           moris::Cell< MSI::Dof_Type > aDofType )
        {
            // get pdofs values for the element
            this->compute_my_pdof_values();

            // get a specific dof type pdofs values
            Matrix< DDRMat > tPdofValues;

            moris::Cell< Matrix< DDRMat > > tPdofValues_Original;

            this->get_my_pdof_values( aDofType, tPdofValues_Original );

            // reshape tCoeffs into the order the cluster expects them
            this->reshape_pdof_values( tPdofValues_Original, tPdofValues );

            // select the required nodal value
            Matrix< IndexMat > tElemVerticesIndices = mMasterInterpolationCell->get_vertex_inds();
            uint tElemNumOfVertices = mMasterInterpolationCell->get_number_of_vertices();

            moris_index tVertexIndex = MORIS_INDEX_MAX;
            for( uint i = 0; i < tElemNumOfVertices; i++ )
            {
                if ( tElemVerticesIndices( i ) == aVertexIndex )
                {
                    tVertexIndex =  i ;
                    break;
                }
            }
            return tPdofValues( tVertexIndex );
        }

//------------------------------------------------------------------------------
    protected:
//------------------------------------------------------------------------------
        /**
         * compute the cluster volume
         */
        real compute_volume();

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
