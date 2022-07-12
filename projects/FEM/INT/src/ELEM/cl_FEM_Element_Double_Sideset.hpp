/*
 * cl_FEM_Element_Double_Sideset.hpp
 *
 *  Created on: May 13, 2019
 *  Author: noel
 */

#ifndef SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_
#define SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_

#include "assert.h"
#include "cl_FEM_Element.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
        class Set;
        //------------------------------------------------------------------------------
        /**
         * \brief Element_Double_Sideset class
         */
        class Element_Double_Sideset : public Element
        {

                //------------------------------------------------------------------------------
            public:
                //------------------------------------------------------------------------------
                /**
                 * constructor
                 * @param[ in ] aLeftIGCell         pointer to mesh cell for master element
                 * @param[ in ] aRightIGCell        pointer to mesh cell for slave element
                 * @param[ in ] aSet                pointer to FEM set to which the elements belong
                 * @param[ in ] aCluster            pointer to FEM cluster to which the elements belong
                 * @param[ in ] aCellIndexInCluster index of the element in the cluster
                 */
                Element_Double_Sideset(
                        mtk::Cell const    * aLeftIGCell,
                        mtk::Cell const    * aRightIGCell,
                        Set                * aSet,
                        Cluster            * aCluster,
                        moris::moris_index   aCellIndexInCluster );

                //------------------------------------------------------------------------------
                /**
                 * destructor
                 */
                ~Element_Double_Sideset();

                //------------------------------------------------------------------------------
                /**
                 * compute residual over the element
                 */
                void compute_residual();

                //------------------------------------------------------------------------------
                /**
                 * compute jacobian over the element
                 */
                void compute_jacobian();

                //------------------------------------------------------------------------------
                /**
                 * compute jacobian and residual over the element
                 */
                void compute_jacobian_and_residual();

                //------------------------------------------------------------------------------
                /**
                 * compute dRdp
                 */
                void compute_dRdp();

                //------------------------------------------------------------------------------
                /**
                 * compute QI
                 */
                void compute_QI();

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdu
                 */
                void compute_dQIdu()
                {
                    MORIS_ERROR( false, "Element_Double_Sideset::compute_dQIdu - not implemented.");
                }

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp
                 */
                void compute_dQIdp_explicit()
                {
                    MORIS_ERROR( false, "Element_Double_Sideset::compute_dQIdp_explicit - not implemented.");
                }

                //------------------------------------------------------------------------------
                /**
                 * compute dRdp and dQIdp
                 */
                void compute_dRdp_and_dQIdp();

                //------------------------------------------------------------------------------
                /**
                 * compute volume over the element
                 * @param[ in ] aIsMaster enum for master or slave
                 */
                real compute_volume( mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

                //------------------------------------------------------------------------------
            protected:

                //------------------------------------------------------------------------------
                /**
                 * initialize the geometry interpolator for the IG master and slave element
                 * @param[ in ] aMasterSideOrdinal side ordinal for the master element
                 * @param[ in ] aSlaveSideOrdinal  side ordinal for the slave element
                 */
                void init_ig_geometry_interpolator(
                        uint aMasterSideOrdinal,
                        uint aSlaveSideOrdinal );

                //------------------------------------------------------------------------------
                /**
                 * initialize the geometry interpolator for the IG master and slave element
                 * @param[ in ] aMasterSideOrdinal side ordinal for the master element
                 * @param[ in ] aSlaveSideOrdinal  side ordinal for the slave element
                 * @param[ in ] aGeoLocalAssembly  matrix with pdv local assembly indices
                 *                                 for master element
                 *                                 ( NumVertexIndices x NumPdvTypes )
                 */
                void init_ig_geometry_interpolator(
                        uint               aMasterSideOrdinal,
                        uint               aSlaveSideOrdinal,
                        Matrix< DDSMat > & aGeoLocalAssembly );
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ELEMENT_DOUBLE_SIDESET_HPP_ */
