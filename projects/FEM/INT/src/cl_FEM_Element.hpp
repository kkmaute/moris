/*
 * cl_FEM_Element.hpp
 *
 *  Created on: Apr 20, 2019
 *      Author: Schmidt
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_HPP_

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

#include "cl_FEM_Set.hpp"   //FEM/INT/src
#include "cl_FEM_Cluster.hpp"   //FEM/INT/src

namespace moris
{
    namespace fem
    {
    class Set;
//------------------------------------------------------------------------------
    /**
     * \brief element class that communicates with the mesh interface
     */
    class Element
    {

    protected:

        //! pointer to master and slave integration cells on mesh
        const mtk::Cell * mMasterCell;
        const mtk::Cell * mSlaveCell;

        moris::moris_index mCellIndexInCluster;

        //! node indices of this element
        //  @node: MTK interface returns copy of vertices. T
        //         storing the indices in private matrix is faster,
        //         but might need more memory
        moris::Matrix< IndexMat > mNodeIndices;

        Set      * mSet     = nullptr;
        Cluster  * mCluster = nullptr;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------
        /**
         * trivial constructor
         */
        Element(){};

//------------------------------------------------------------------------------
        /**
         * constructor for master only
         * @param[ in ] aCell               a mesh cell pointer
         * @param[ in ] aSet                a fem set pointer
         * @param[ in ] aCluster            a mesh cluster pointer
         * @param[ in ] aCellIndexInCluster the index of the cell within the cluster
         */
        Element( const mtk::Cell    * aCell,
                 Set                * aSet,
                 Cluster            * aCluster,
                 moris::moris_index   aCellIndexInCluster ) : mSet( aSet ),
                                                              mCluster( aCluster )
        {
            // fill the cell index in cluster
            mCellIndexInCluster = aCellIndexInCluster;

            // fill the bulk mtk::Cell pointer //FIXME
            mMasterCell = aCell;
        };

//------------------------------------------------------------------------------
        /**
         * constructor for slave only
         * @param[ in ] aMasterCell         a master mesh cell pointer
         * @param[ in ] aSlaveCell          a slave mesh cell pointer
         * @param[ in ] aSet                a fem set pointer
         * @param[ in ] aCluster            a mesh cluster pointer
         * @param[ in ] aCellIndexInCluster the index of the cell within the cluster
         */
        Element( const mtk::Cell    * aMasterCell,
                 const mtk::Cell    * aSlaveCell,
                 Set                * aSet,
                 Cluster            * aCluster,
                 moris::moris_index   aCellIndexInCluster ) : mSet( aSet ),
                                                              mCluster( aCluster )
        {
            // fill the cell index in cluster
            mCellIndexInCluster = aCellIndexInCluster;

            // fill the master and slave cell pointers
            mMasterCell = aMasterCell;
            mSlaveCell  = aSlaveCell;
        };

//------------------------------------------------------------------------------
        /**
         * trivial destructor
         */
        virtual ~Element(){};

//------------------------------------------------------------------------------
        /**
         * compute jacobian
         */
        virtual void compute_jacobian() = 0;

//------------------------------------------------------------------------------
        /**
         * compute residual
         */
        virtual void compute_residual() = 0;

//------------------------------------------------------------------------------
        /**
         * compute jacobian and residual
         */
        virtual void compute_jacobian_and_residual() = 0;

//------------------------------------------------------------------------------
        /**
         * compute quantity of interest
         * @param[ in ] aQIComputeType an enum for computation type
         */
        void compute_quantity_of_interest( fem::QI_Compute_Type aQIComputeType )
        {
            switch ( aQIComputeType )
            {
                case ( fem::QI_Compute_Type::GLOBAL ) :
                {
                    this->compute_quantity_of_interest_global();
                    break;
                }
                case ( fem::QI_Compute_Type::NODAL ) :
                {
                    this->compute_quantity_of_interest_nodal();
                    break;
                }
                case ( fem::QI_Compute_Type::ELEMENTAL ) :
                {
                    this->compute_quantity_of_interest_elemental();
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "Element::compute_quatity_of_interest - unknow QI compute type." );
                    break;
                }
            }
        }

//------------------------------------------------------------------------------
        /**
         * compute quantity of interest in a global way
         */
        void compute_quantity_of_interest_global()
        {
            MORIS_ERROR( false, "Element::compute_quantity_of_interest_global - this function does nothing." );
        }

//------------------------------------------------------------------------------
        /**
         * compute quantity of interest in a nodal way
         */
        void compute_quantity_of_interest_nodal()
        {
            MORIS_ERROR( false, "Element::compute_quantity_of_interest_nodal - this function does nothing." );
        }

//------------------------------------------------------------------------------
        /**
         * compute quantity of interest in an elemental way
         */
        void compute_quantity_of_interest_elemental()
        {
            MORIS_ERROR( false, "Element::compute_quantity_of_interest_elemental - this function does nothing." );
        }

//------------------------------------------------------------------------------
        /**
         * compute volume of the integration element
         * @param[ in ] aIsMaster enum master or slave
         */
        real compute_volume( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
        {
            //get number of integration points
            uint tNumOfIntegPoints = mSet->get_number_of_integration_points();

            // init volume
            real tVolume = 0;

            // loop over integration points
            for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // set integration point for geometry interpolator
                mSet->get_IG_geometry_interpolator()->set_space_time( mSet->get_integration_points().get_column( iGP ) );

                // compute and add integration point contribution to volume
                tVolume += mSet->get_IG_geometry_interpolator()->det_J() * mSet->get_integration_weights()( iGP );
            }

            // return the volume value
            return tVolume;
        }

//------------------------------------------------------------------------------
        /**
          * compute size of the integration element
          * @param[ in ] aSpaceDim space dimension
          * FIXME should be held by the FEM Set
          * @param[ in ] aIsMaster enum master or slave
          */
         real compute_size( uint              aSpaceDim,
                            mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER )
         {
             //evaluate the volume
             real tVolume = this->compute_volume( aIsMaster );

             // compute the element size
             switch ( aSpaceDim )
             {
                 case ( 3 ):
                     return std::pow( 6 * tVolume / M_PI, 1.0 / 3.0 );

                 case ( 2 ):
                     return std::pow( 4 * tVolume / M_PI, 1.0 / 2.0 );

                 case ( 1 ):
                     return tVolume;

                 default:
                     MORIS_ERROR( false, "Element::compute_size - space dimension can only be 1, 2, or 3. ");
                     return tVolume;
             }
         }

//------------------------------------------------------------------------------
    };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ELEMENT_HPP_ */
