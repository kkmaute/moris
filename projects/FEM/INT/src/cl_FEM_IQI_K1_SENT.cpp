/*
 * cl_FEM_IQI_K1_SENT.cpp
 *
 *  Created on: Mar 13, 2020
 *      Author: sonne
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_K1_SENT.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
            IQI_K1_SENT::IQI_K1_SENT()
            {
                // set IQI type
                mIQIType = vis::Output_Type::K1_SENT;

                // set fem IQI type
                mFEMIQIType = fem::IQI_Type::K1_SENT;

                // set size for the constitutive model pointer cell
                mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

                // populate the constitutive map
                mConstitutiveMap[ "ElastLinIso" ] = IQI_Constitutive_Type::ELAST_LIN_ISO;

                // set size for the property pointer cell
                mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

                // populate the property map
                mPropertyMap[ "CRACK_LENGTH" ] = IQI_Property_Type::CRACK_LENGTH;
                mPropertyMap[ "NUM_NODES" ]    = IQI_Property_Type::NUM_NODES;
                mPropertyMap[ "E" ]            = IQI_Property_Type::E;
                mPropertyMap[ "NU" ]           = IQI_Property_Type::NU;
                mPropertyMap[ "PLAIN_TYPE" ]   = IQI_Property_Type::PLAIN_TYPE;

            }
//------------------------------------------------------------------------------
            void IQI_K1_SENT::compute_QI( Matrix< DDRMat > & aQI )
            {
                /*
                 * TODO:
                 *      - loop over both the top and the bottom nodes
                 *      - get plain type from CM? or input as a parameter (plain stress/strain)
                 */

                real tPi = 3.14159265358979323846;

                /* -------------------- get node data -------------------- */
                // get node displacement values
                uint tMeshIndex = 0;

                Matrix< DDRMat > tElemVals;     // dummy matrix
                Matrix< DDRMat > tNodeVals;     // node displacements
                real             tGlobVal;      // dummy global value
                mSet->compute_quantity_of_interest( tMeshIndex, &tElemVals, &tNodeVals, &tGlobVal, vis::Output_Type::DOF, vis::Field_Type::NODAL );

print(tNodeVals,"tNodeVals");
                // get node coordinates
                Cell< Matrix< DDRMat > > tAllCoords;                                                        // cell of all node coordinate vectors
                moris::Cell< mtk::Cluster const * > tTempClusters = mSet->get_clusters_on_set();

                uint tNumClusters = tTempClusters.size();

                for(uint iClust=0; iClust<tNumClusters; iClust++)                                           // loop over all clusters on set
                {
                    moris::mtk::Cell const & tIPCell = tTempClusters(iClust)->get_interpolation_cell();

                    moris::Cell< moris::mtk::Vertex * > tVertices = tIPCell.get_vertex_pointers();
                    uint tNumVerts = tVertices.size();

                    uint tOldSize = tAllCoords.size();

                    tAllCoords.resize( tOldSize + tNumVerts );

                    for(uint iVert=0; iVert<tNumVerts; iVert++)                                             // loop over all vertices on the cluster
                    {
                        tAllCoords( tOldSize + iVert ) = tVertices( iVert )->get_coords();                  // fill the element with coordinate values
                    }
                }
print(tAllCoords,"tAllCoords");
                /* -------------------- end get node data ---------------- */
                uint tCrackLengthIndex = static_cast< uint >( IQI_Property_Type::CRACK_LENGTH );                // get property index for crack length
                real tCrackLength      = mMasterProp( tCrackLengthIndex )->val()(0,0);

                uint tNumNodesIndex = static_cast< uint >( IQI_Property_Type::NUM_NODES );                      // get property index for number of nodes to integrate over
                real tNumNodes      = mMasterProp( tNumNodesIndex )->val()(0,0);

                Matrix< DDRMat > tAMat(3,3, 0.0);
                Matrix< DDRMat > tBVec(3,1, 0.0);

                for( uint ii=0; ii<1; ii++ )        // loop over top and bottom nodes ( needs to be implemented for bottom nodes )
                {
                    for( uint i=2; i<(uint)tNumNodes; i++ )     // loop over specified number of nodes along crack line ( the number 100 should be pulled from a parameter )
                    {
                        if( ii==0 )
                        {
                            /*  for the top nodes   */
                            real tR  = tCrackLength - tAllCoords( i )(0);                // radius defined as distance in x-direction from the crack tip
                            real tDR = tAllCoords( i )(0) - tAllCoords( i-1 )(0);
                            real tUx = tNodeVals( i,0 );    // TODO: make sure this is how the information is being returned in the compute_IQI function above

                            tAMat(0,0) += (tR/2/tPi)*tDR;
                            tAMat(0,1) += std::sqrt(tR/2/tPi)*tDR;
                            tAMat(0,2) += std::sqrt(tR/2/tPi)*tR*tDR;

                            tAMat(1,0) += std::sqrt(tR/2/tPi)*tDR;
                            tAMat(1,1) += tDR;
                            tAMat(1,2) += tR*tDR;

                            tAMat(2,0) += std::sqrt(tR/2/tPi)*tR*tDR;
                            tAMat(2,1) += tR*tDR;
                            tAMat(2,2) += (tR*tR)*tDR;

                            tBVec(0) += std::sqrt(tR/2/tPi)*tUx*tDR;
                            tBVec(1) += tUx*tDR;
                            tBVec(2) += tUx*tR*tDR;
                        }
                        if( ii==1 )
                        {
                            /*  for the bottom nodes   */
                        }
                    }
                }

                /* solve Ax=b problem to get coefficients and use them to compute K1 */
                Matrix< DDRMat > tCVec = inv(tAMat)*tBVec;

                /*
                 * case: plain stress
                 *          Eprime = E
                 * case: plain strain
                 *          Eprime = E/(1-nu^2)
                 */
                uint tEIndex = static_cast< uint >( IQI_Property_Type::E );                 // get property index for Young's Modulus
                real tEprime = mMasterProp( tEIndex )->val()(0,0);                          // FIXME: get plain type from CM? currently assuming plain stress....

                aQI.resize(1,1);
                /* compute K1 from displacement correlation */
                aQI(0,0) = tEprime*tCVec(0)/4;
print(aQI,"K1");
            }
//------------------------------------------------------------------------------
            void IQI_K1_SENT::compute_dQIdu(  Matrix< DDRMat > & adQIdDof  )
            {
                MORIS_ERROR( false, "IQI_K1_SENT::compute_dQIdu() - this function does nothing for K1_SENT IQI type" );
            }
    }   // end fem namespace
}       // end moris namespace
