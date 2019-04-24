/*
 * cl_Equation_Object.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_OBJECT_HPP_
#define SRC_FEM_CL_EQUATION_OBJECT_HPP_

#include <memory>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "fn_trans.hpp"
#include "op_times.hpp"

#include "cl_MSI_Pdof_Host.hpp"
namespace moris
{
class Dist_Vector;
    namespace fem
    {
        class Node_Base;
        class Element;
    }
    namespace MSI
    {
        class Pdof;
        class Pdof_Host;
        class Equation_Object
        {
//-------------------------------------------------------------------------------------------------
        protected:
//-------------------------------------------------------------------------------------------------
            moris::Cell< fem::Node_Base * >         mNodeObj;
            moris::Cell< Pdof_Host * >              mMyPdofHosts;       // Pointer to the pdof hosts of this equation object

            moris::Cell< enum Dof_Type >            mEqnObjDofTypeList; // List of dof types of this equation obj
            moris::Cell< Pdof* >                    mFreePdofs;         // List of the pdof pointers of this equation obj

            Matrix< DDSMat >                        mUniqueAdofList;    // Unique adof list for this equation object
            moris::map < moris::uint, moris::uint > mUniqueAdofMap;     // Map to

            //! weak BCs of element
            Matrix< DDRMat > mNodalWeakBCs;

            Matrix< DDRMat > mResidual;
            Matrix< DDRMat > mJacobian;

            Matrix< DDRMat > mPdofValues;

            Dist_Vector * mSolVec = nullptr;

            Model_Solver_Interface * mModelSolverInterface = nullptr;

            moris::uint mEqnObjInd;

            // sideset information //FIXME Side ordinals are not part of the equation object
            Matrix< IndexMat > mListOfSideOrdinals;
            Matrix< IndexMat > mListOfTimeOrdinals;

            Matrix< DDRMat >mTime;

            friend class fem::Element;

//-------------------------------------------------------------------------------------------------
        public:
//-------------------------------------------------------------------------------------------------

            Equation_Object() {};

//-------------------------------------------------------------------------------------------------
            Equation_Object( const moris::Cell< fem::Node_Base * > & aNodeObjs );

//-------------------------------------------------------------------------------------------------

            virtual ~Equation_Object(){};

//-------------------------------------------------------------------------------------------------

            void set_time( const Matrix< DDRMat > & aTime )
            {
                mTime = aTime;
            }

//-------------------------------------------------------------------------------------------------

            void set_model_solver_interface_pointer( Model_Solver_Interface * aModelSolverInterface)
            {
                mModelSolverInterface = aModelSolverInterface;
            };

//-------------------------------------------------------------------------------------------------

            Matrix< DDRMat > & get_pdof_values( )
            {
                this->get_my_pdof_values();

                return mPdofValues;
            };

//-------------------------------------------------------------------------------------------------
            /**
             * @brief Get function to get the dof types used by this equation object. This function is tested by the test [Dof_Mgn_create_unique_dof_type_list]
             * [Dof_Mgn_create_unique_dof_type_map_matrix]
             *
             * @param[in] aDofType   List of dof types.
             *
             */
            void get_dof_types( moris::Cell< enum Dof_Type > & aDofType )
            {
                aDofType = mEqnObjDofTypeList;

//                // get the number of groups of dof types
//                uint tNumOfDofTypeGropups = mEqnObjDofTypeList.size();
//
//                // get the total number of dof types
//                uint tCounter = 0;
//                for ( uint i = 0; i < tNumOfDofTypeGropups; i++ )
//                {
//                    tCounter = tCounter + mEqnObjDofTypeList( i ).size();
//                }
//
//                // set the size of aDofType
//                aDofType.resize( tCounter );
//
//                // loop over the groups of dof types
//                tCounter = 0;
//                for ( uint i = 0; i < tNumOfDofTypeGropups; i++ )
//                {
//                    Cell< MSI::Dof_Type > tDofTypeGroup = mEqnObjDofTypeList( i );
//                    uint tNumOfDofTypesInGroupI = tDofTypeGroup.size();
//
//                    for( uint j = 0; j < tNumOfDofTypesInGroupI; j++ )
//                    {
//                        aDofType( tCounter ) = tDofTypeGroup( j );
//                        tCounter++;
//                    }
//                }
            }
//-------------------------------------------------------------------------------------------------
            /**
             * @brief Returns the number of nodes, elements and ghosts related to this equation object.
             *
             */
            // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
            moris::uint get_num_pdof_hosts() { return mNodeObj.size(); }

//------------------------------------------------------------------------------------------------
            /**
             * @brief Returns the maximal pdof host (node) index of this equation object
             *
             */
            moris::uint get_max_pdof_hosts_ind();

//-------------------------------------------------------------------------------------------------
            /**
             * @brief Creates the pdof hosts of this equation object, if not created earlier, and puts them into the local pdof host list. This function is tested by the test [Eqn_Obj_create_pdof_host]
             *
             * @param[in] aNumUsedDofTypes   Number of globally used dof types
             * @param[in] aPdofTypeMap       Map which maps the dof type enum values to a consecutive list of dof type indices.
             * @param[in] aPdofHostList      List of pdof hosts.
             *
             */
            void create_my_pdof_hosts( const moris::uint                  aNumUsedDofTypes,
                                       const Matrix< DDSMat >           & aPdofTypeMap,
                                       const Matrix< DDUMat >           & aTimePerDofType,
                                             moris::Cell< Pdof_Host * > & aPdofHostList );

//-------------------------------------------------------------------------------------------------
            /**
             * @brief This function creates a list of pdof pointers related to this equation object. This function is tested by the test [Eqn_Obj_create_my_pdof_list]
             * [Dof_Mgn_create_unique_dof_type_map_matrix]
             *
             */
            void create_my_pdof_list();

//-------------------------------------------------------------------------------------------------
            /**
             * @brief This function creates a unique list of adofs Ids corresponding to this equation object. This function is tested by the test [Eqn_Obj_create_my_list_of_adof_ids]
             *
             */
            void create_my_list_of_adof_ids();

//-------------------------------------------------------------------------------------------------
            /**
             * @brief This function creates a map relating the adof ids to the positions for this equation object . This function is tested by the test [Eqn_Obj_create_adof_map]
             *
             */
            void set_unique_adof_map();

//-------------------------------------------------------------------------------------------------
            /**
             * @brief This function creates a PADofMap witch can be used to for a calculation from pdofs to adofs . This function is tested by the test [Eqn_Obj_PADofMap]
             *
             */
            void build_PADofMap( Matrix< DDRMat > & aPADofMap );

//-------------------------------------------------------------------------------------------------

            /**
             * @brief Get function for the pdof values of this particular equation object
             *
             */
            void get_my_pdof_values( );

//-------------------------------------------------------------------------------------------------

            /**
             * @brief Get function for the pdof values of this particular equation object. get_my_pdof_values() has to be called first to initialize.
             *
             * @param[in] aRequestedDofTypes      List of requested dof types
             * @param[in] aRequestedPdofValues    Reference to the matrix of requested pdof values
             */
            void get_my_pdof_values( const moris::Cell< enum Dof_Type > & aRequestedDofTypes,
                                           Matrix< DDRMat >             & aRequestedPdofValues);

//-------------------------------------------------------------------------------------------------

            void set_vector_entry_number_of_pdof();

//-------------------------------------------------------------------------------------------------

            void get_egn_obj_jacobian( Matrix< DDRMat > & aEqnObjMatrix,
                                       Dist_Vector      * aSolutionVector )
            {
                mSolVec = aSolutionVector;

                Matrix< DDRMat > tTMatrix;
                this->build_PADofMap( tTMatrix );

                this->compute_jacobian();

//                print( tTMatrix, "tTMatrix" );
//                print( mJacobian, "mJacobian" );

                aEqnObjMatrix = trans( tTMatrix ) * mJacobian * tTMatrix ;

            };

//-------------------------------------------------------------------------------------------------

            void get_equation_obj_residual( Matrix< DDRMat > & aEqnObjRHS,
                                            Dist_Vector * aSolutionVector );

//-------------------------------------------------------------------------------------------------
            void get_equation_obj_dof_ids( Matrix< DDSMat > & aEqnObjAdofId )
            {
                aEqnObjAdofId = mUniqueAdofList;
            };
//-------------------------------------------------------------------------------------------------
            /**
             * returns a moris::Mat with indices of vertices that are connected to this element
             */
            moris_index get_node_index( const moris_index aElementLocalNodeIndex ) const ;

//-------------------------------------------------------------------------------------------------

            virtual Matrix< DDSMat > get_adof_indices()
            {
                MORIS_ERROR( false, "this function does nothing");
                return Matrix< DDSMat >(0,0);
            }

//-------------------------------------------------------------------------------------------------

            virtual void compute_jacobian()
            {
                MORIS_ERROR( false, "this function does nothing");
            }

//-------------------------------------------------------------------------------------------------
            virtual void compute_residual()
            {
                MORIS_ERROR( false, "this function does nothing");
            }

//-------------------------------------------------------------------------------------------------
            virtual void compute_jacobian_and_residual()
            {
                MORIS_ERROR( false, "this function does nothing");
            }

//-------------------------------------------------------------------------------------------------
            virtual moris::real compute_integration_error( moris::real (*aFunction)( const Matrix< DDRMat > & aPoint ) )
            {
                MORIS_ERROR( false, "this function does nothing");
                return 0.0;
            }

//-------------------------------------------------------------------------------------------------
            virtual moris::real compute_element_average_of_scalar_field()
            {
                MORIS_ERROR( false, "this function does nothing");
                return 0.0;
            }

//-------------------------------------------------------------------------------------------------
            /**
             * return Neumann boundary conditions, writable version
             */
            virtual Matrix< DDRMat > & get_weak_bcs()
            {
                return mNodalWeakBCs;
            }

//-------------------------------------------------------------------------------------------------
            /**
             * return Neumann boundary conditions, const version
             */
            const Matrix< DDRMat > & get_weak_bcs() const
            {
                return mNodalWeakBCs;
            }

//-------------------------------------------------------------------------------------------------
            /**
             * set the list of side ordinals
             */
            void set_list_of_side_ordinals( const Matrix< IndexMat > & aListOfSideOrdinals )
            {
                mListOfSideOrdinals = aListOfSideOrdinals;
            }

//-------------------------------------------------------------------------------------------------
            /**
             * set the list of time ordinals
             */
            void set_list_of_time_ordinals( const Matrix< IndexMat > & aListOfTimeOrdinals )
            {
                mListOfTimeOrdinals = aListOfTimeOrdinals;
            }

//-------------------------------------------------------------------------------------------------
            /**
             * how many nodes are connected to this element
             */
            uint get_num_nodes() const
            {
                return mNodeObj.size();
            }

//-------------------------------------------------------------------------------------------------

            virtual moris::real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                                                       moris::Cell< MSI::Dof_Type > aDofType )
            {
                MORIS_ERROR( false, "Equation_Object::get_element_nodal_pdof_value - this function does nothing");
                return 0.0;
            }
        };
    }
}

#endif /* SRC_FEM_CL_EQUATION_OBJECT_HPP_ */
