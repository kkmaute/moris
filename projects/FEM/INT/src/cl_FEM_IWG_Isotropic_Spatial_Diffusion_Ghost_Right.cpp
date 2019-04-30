
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost_Right.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Ghost_Right::IWG_Isotropic_Spatial_Diffusion_Ghost_Right()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };

            // FIXME set a penalty
            mGammaGhost = 0.0001;

            // FIXME set mesh parameter
            // Mesh parameter or parameters which it depends on must be fed to IWG from outside
            mMeshParameter = 1;

            // FIXME set order of Shape functions
            // Order must be fed to IWG from outside
            mOrder = 2;

            //FIXME forced diffusion parameter
            //      forced dimensions for 3D
            eye( mSpaceDim, mSpaceDim, mKappa );
            mKappa = 1.0 * mKappa;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::compute_residual
            ( Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check, if order of shape functions is supported
        	MORIS_ERROR( mOrder < 3, " Ghost stabilization for this order of shape fncts. not supported yet. " );

        	// set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // compute the residual r_T
            uint tResSize = tTemp->get_number_of_space_time_coefficients();
            aResidual.set_size( tResSize, 1, 0.0);

            // --------- Ghost for linear shape functions ---------
            if (mOrder == 1) {
                aResidual =   aResidual
            		        - mGammaGhost * mMeshParameter // * mKappa                            // scaling parameters
            		        * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->gradx( 1 );  // matrices & vectors
            }

            // --------- Ghost for quadratic shape functions ---------
            else if (mOrder == 2) {

            	Matrix <DDRMat> grad2x = tTemp->gradx(2);
            	Matrix <DDRMat> B2x = tTemp->eval_d2Ndx2();
            	uint numOfShapeFunctions = B2x.n_cols();

            	Matrix <DDRMat> B2xNormal( 3, numOfShapeFunctions, 0.0 );
            	for(uint i = 0; i < numOfShapeFunctions; i++){
            		B2xNormal(0,i) = B2x(0,i) * mNormal(0) + B2x(3,i) * mNormal(1) + B2x(5,i) * mNormal(2);
            		B2xNormal(1,i) = B2x(3,i) * mNormal(0) + B2x(1,i) * mNormal(1) + B2x(4,i) * mNormal(2);
            		B2xNormal(2,i) = B2x(5,i) * mNormal(0) + B2x(4,i) * mNormal(1) + B2x(2,i) * mNormal(2);
            	}

            	Matrix <DDRMat> grad2xNormal( 3, 1, 0.0 );
            	grad2xNormal(0) = grad2x(0) * mNormal(0) + grad2x(3) * mNormal(1) + grad2x(5) * mNormal(2);
        		grad2xNormal(1) = grad2x(3) * mNormal(0) + grad2x(1) * mNormal(1) + grad2x(4) * mNormal(2);
        		grad2xNormal(2) = grad2x(5) * mNormal(0) + grad2x(4) * mNormal(1) + grad2x(2) * mNormal(2);

            	aResidual =   aResidual
            	            - mGammaGhost * mMeshParameter // * mKappa                            // scaling parameters
            	            * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->gradx( 1 )   // matrices & vectors
							- mGammaGhost * (std::pow(mMeshParameter,3)) // * mKappa              // scaling parameters (p=2)
							* trans(B2xNormal) * grad2xNormal;                                    // matrices & vectors (p=2)
            }

            // --------- Ghost for other shape functions not supported (yet) ---------
            else {
            	MORIS_ERROR( false, " Ghost stabilization for this order of shape fncts. not supported yet. " );
            }

        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::compute_jacobian
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check, if order of shape functions is supported
        	MORIS_ERROR( mOrder < 3, " Ghost stabilization for this order of shape fncts. not supported yet. " );

        	// set field interpolator
            Field_Interpolator* tTemp  = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            uint tNumOfBases = tTemp->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );

            // --------- Ghost for linear shape functions ---------
            if (mOrder == 1) {
            	 aJacobians( 0 ) =    aJacobians( 0 )
            		        - mGammaGhost * mMeshParameter // * mKappa                      // scaling parameters
            		        * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->Bx();  // matrices & vectors
            }

            // --------- Ghost for quadratic shape functions ---------
            else if (mOrder == 2) {

            	Matrix <DDRMat> B2x = tTemp->eval_d2Ndx2();
            	uint numOfShapeFunctions = B2x.n_cols();

            	Matrix <DDRMat> B2xNormal( 3, numOfShapeFunctions, 0.0 );
            	for(uint i = 0; i < numOfShapeFunctions; i++){
            		B2xNormal(0,i) = B2x(0,i) * mNormal(0) + B2x(3,i) * mNormal(1) + B2x(5,i) * mNormal(2);
            		B2xNormal(1,i) = B2x(3,i) * mNormal(0) + B2x(1,i) * mNormal(1) + B2x(4,i) * mNormal(2);
            		B2xNormal(2,i) = B2x(5,i) * mNormal(0) + B2x(4,i) * mNormal(1) + B2x(2,i) * mNormal(2);
            	}

            	aJacobians( 0 ) =   aJacobians( 0 )
            	                  - mGammaGhost * mMeshParameter // * mKappa                    // scaling parameters
            	                  * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->Bx() // matrices & vectors
							      - mGammaGhost * (std::pow(mMeshParameter,3)) // * mKappa      // scaling parameters (p=2)
							      * trans(B2xNormal) * B2xNormal;                               // matrices & vectors (p=2)
            }

            // --------- Ghost for other shape functions not supported (yet) ---------
            else {
            	MORIS_ERROR( false, " Ghost stabilization for this order of shape fncts. not supported yet. " );
            }

        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::compute_jacobian_and_residual
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // check, if order of shape functions is supported
        	MORIS_ERROR( mOrder < 3, " Ghost stabilization for this order of shape fncts. not supported yet. " );

        	// set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            //______________________________________________//
            //////////////////// RESIDUAL ////////////////////
            //----------------------------------------------//

            // compute the residual r_T
            uint tResSize = tTemp->get_number_of_space_time_coefficients();
            aResidual.set_size( tResSize, 1, 0.0);

            // --------- Ghost for linear shape functions ---------
            if (mOrder == 1) {
                aResidual =   aResidual
            		        - mGammaGhost * mMeshParameter // * mKappa                            // scaling parameters
            		        * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->gradx( 1 );  // matrices & vectors
            }

            // --------- Ghost for quadratic shape functions ---------
            else if (mOrder == 2) {

            	Matrix <DDRMat> grad2x = tTemp->gradx(2);
            	Matrix <DDRMat> B2x = tTemp->eval_d2Ndx2();
            	uint numOfShapeFunctions = B2x.n_cols();

            	Matrix <DDRMat> B2xNormal( 3, numOfShapeFunctions, 0.0 );
            	for(uint i = 0; i < numOfShapeFunctions; i++){
            		B2xNormal(0,i) = B2x(0,i) * mNormal(0) + B2x(3,i) * mNormal(1) + B2x(5,i) * mNormal(2);
            		B2xNormal(1,i) = B2x(3,i) * mNormal(0) + B2x(1,i) * mNormal(1) + B2x(4,i) * mNormal(2);
            		B2xNormal(2,i) = B2x(5,i) * mNormal(0) + B2x(4,i) * mNormal(1) + B2x(2,i) * mNormal(2);
            	}

            	Matrix <DDRMat> grad2xNormal( 3, 1, 0.0 );
            	grad2xNormal(0) = grad2x(0) * mNormal(0) + grad2x(3) * mNormal(1) + grad2x(5) * mNormal(2);
        		grad2xNormal(1) = grad2x(3) * mNormal(0) + grad2x(1) * mNormal(1) + grad2x(4) * mNormal(2);
        		grad2xNormal(2) = grad2x(5) * mNormal(0) + grad2x(4) * mNormal(1) + grad2x(2) * mNormal(2);

            	aResidual =   aResidual
            	            - mGammaGhost * mMeshParameter // * mKappa                            // scaling parameters
            	            * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->gradx( 1 )   // matrices & vectors
							- mGammaGhost * (std::pow(mMeshParameter,3)) // * mKappa              // scaling parameters (p=2)
							* trans(B2xNormal) * grad2xNormal;                                    // matrices & vectors (p=2)
            }

            // --------- Ghost for other shape functions not supported (yet) ---------
            else {
            	MORIS_ERROR( false, " Ghost stabilization for this order of shape fncts. not supported yet. " );
            }

            //______________________________________________//
            //////////////////// JACOBIAN ////////////////////
            //----------------------------------------------//

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            uint tNumOfBases = tTemp->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );

            // --------- Ghost for linear shape functions ---------
            if (mOrder == 1) {
            	 aJacobians( 0 ) =    aJacobians( 0 )
            		        - mGammaGhost * mMeshParameter // * mKappa                      // scaling parameters
            		        * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->Bx();  // matrices & vectors
            }

            // --------- Ghost for quadratic shape functions ---------
            else if (mOrder == 2) {

            	Matrix <DDRMat> B2x = tTemp->eval_d2Ndx2();
            	uint numOfShapeFunctions = B2x.n_cols();

            	Matrix <DDRMat> B2xNormal( 3, numOfShapeFunctions, 0.0 );
            	for(uint i = 0; i < numOfShapeFunctions; i++){
            		B2xNormal(0,i) = B2x(0,i) * mNormal(0) + B2x(3,i) * mNormal(1) + B2x(5,i) * mNormal(2);
            		B2xNormal(1,i) = B2x(3,i) * mNormal(0) + B2x(1,i) * mNormal(1) + B2x(4,i) * mNormal(2);
            		B2xNormal(2,i) = B2x(5,i) * mNormal(0) + B2x(4,i) * mNormal(1) + B2x(2,i) * mNormal(2);
            	}

            	aJacobians( 0 ) =   aJacobians( 0 )
            	                  - mGammaGhost * mMeshParameter // * mKappa                    // scaling parameters
            	                  * trans(tTemp->Bx()) * mNormal * trans(mNormal) * tTemp->Bx() // matrices & vectors
							      - mGammaGhost * (std::pow(mMeshParameter,3)) // * mKappa      // scaling parameters (p=2)
							      * trans(B2xNormal) * B2xNormal;                               // matrices & vectors (p=2)                                      // matrices & vectors (p=2)
            }

            // --------- Ghost for other shape functions not supported (yet) ---------
            else {
            	MORIS_ERROR( false, " Ghost stabilization for this order of shape fncts. not supported yet. " );
            }


        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::set_interpolation_order ( uint aOrder )
        {
        	if ( (aOrder == 1) || (aOrder == 2) ){
        		mOrder = aOrder;
        	}

        	else {
        		MORIS_ERROR( false, "Error in method .set_interpolation_order() \nGhost stabilization for this order of shape functions not supported yet. " );
        	}
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Ghost_Right::set_penalty_factor ( real aGamma )
        {
            mGammaGhost = aGamma;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
