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

#include "cl_MTK_Enums.hpp"                 //FEM/INT/src
#include "cl_MTK_Vertex.hpp"      //MTK/src

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
    }
    namespace fem
    {
        class Node_Base;
        class Element;
    }
    namespace fem
    {
        class Cluster;

    }
    namespace vis
    {
        enum class Output_Type;
        enum class Field_Type;
    }
    namespace MSI
    {
        class Pdof;
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
                moris::Cell< moris::Cell< fem::Node_Base * > >     mNodeObj;
                moris::Cell< moris::Cell< Pdof_Host * > >          mMyPdofHosts;       // Pointer to the pdof hosts of this equation object

                moris::Cell< Pdof* >                               mFreePdofs;         // List of the pdof pointers of this equation obj
                moris::Cell< moris::Cell< moris::Cell< Pdof* > > > mFreePdofList;      // FIXME list of free pdofs ordered after their dof type . mFreePdofs or mFreePdofList should be deleted

                Matrix< DDSMat >                                   mUniqueAdofList;    // Unique adof list for this equation object
                moris::Cell< moris::Cell< Matrix< DDSMat > > >     mUniqueAdofTypeList;
                moris::map < moris::uint, moris::uint >            mUniqueAdofMap;     // Map to

                moris::Cell< moris::Cell< moris::map < moris::uint, moris::uint > > > mUniqueAdofMapList;     // Map to

                //! weak BCs of element FIXME
                Matrix< DDRMat > mNodalWeakBCs;

                //! actual pdof values. Cells are for different multi-vectors
                moris::Cell< Matrix< DDRMat > > mPdofValues;

                //! previous pdof values
                moris::Cell< Matrix< DDRMat > > mPreviousPdofValues;

                //! adjoint pdof values
                moris::Cell< Matrix< DDRMat > > mAdjointPdofValues;

                moris::uint mEqnObjInd;

                //            Matrix< DDRMat > mTime;
                //            Matrix< DDRMat > mPrevTime;

                Equation_Set * mEquationSet;

                moris::uint mNumPdofSystems = 0;

                friend class fem::Element;
                friend class fem::Cluster;

                //-------------------------------------------------------------------------------------------------
            public:
                //-------------------------------------------------------------------------------------------------

                Equation_Object(){};

                Equation_Object( Equation_Set * aElementBlock ) : mEquationSet( aElementBlock )
                {};

                //-------------------------------------------------------------------------------------------------
                Equation_Object( const moris::Cell < moris::Cell< fem::Node_Base * > > & aNodeObjs );

                //-------------------------------------------------------------------------------------------------

                virtual ~Equation_Object(){};

                //-------------------------------------------------------------------------------------------------

                void set_time( Matrix< DDRMat > & aTime );

                //-------------------------------------------------------------------------------------------------

                Matrix< DDRMat > & get_time();

                //-------------------------------------------------------------------------------------------------

                Matrix< DDRMat > & get_previous_time();

                //-------------------------------------------------------------------------------------------------
                Cell< Matrix< DDRMat > > & get_pdof_values( )
                {
                    this->compute_my_pdof_values();

                    return mPdofValues;
                };

                //-------------------------------------------------------------------------------------------------
                /**
                 * @brief Returns the number of nodes, elements and ghosts related to this equation object. This function is only for unit test purposes.
                 *
                 */
                // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
                moris::uint get_num_pdof_hosts()
                {
                    moris::uint tNumPdofHosts = 0;
                    for( uint Ik = 0; Ik < mNodeObj.size(); Ik++ )
                    {
                        tNumPdofHosts = tNumPdofHosts + mNodeObj( Ik ).size();
                    }
                    return tNumPdofHosts;
                }

                //------------------------------------------------------------------------------------------------

                /**
                 * @brief Returns the maximal pdof host (node) index of this equation object
                 *
                 */
                moris::uint get_max_pdof_hosts_ind();

                //-------------------------------------------------------------------------------------------------

                /**
                 * @brief Creates the pdof hosts of this equation object, if not created earlier, and puts them into the local pdof host list.
                 *  This function is tested by the test [Eqn_Obj_create_pdof_host]
                 *
                 * @param[in] aNumUsedDofTypes   Number of globally used dof types
                 * @param[in] aPdofTypeMap       Map which maps the dof type enum values to a consecutive list of dof type indices.
                 * @param[in] aPdofHostList      List of pdof hosts.
                 *
                 */
                void create_my_pdof_hosts(
                        const moris::uint                  aNumUsedDofTypes,
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

                void build_PADofMap_list( Cell< Cell< Matrix< DDRMat > > > & aPADofMap );

                void build_PADofMap_1( Matrix< DDRMat > & aPADofMap );

                //-------------------------------------------------------------------------------------------------

                /**
                 * @brief Compute function for the pdof values of this particular equation object
                 *
                 */
                void compute_my_pdof_values( );

                //-------------------------------------------------------------------------------------------------

                /**
                 * @brief Compute function for the previous pdof values of this particular equation object
                 *
                 */
                void compute_previous_pdof_values( );

                //-------------------------------------------------------------------------------------------------

                /**
                 * @brief Compute function for the labda values
                 *
                 */
                void compute_my_adjoint_values( );

                //-------------------------------------------------------------------------------------------------

                /**
                 * @brief Get function for the pdof values of this particular equation object.
                 * get_my_pdof_values() has to be called first to initialize.
                 * @param[ in ] aPdofValues             All pdof values of this equation object
                 * @param[ in ] aRequestedDofTypes      List of requested dof types
                 * @param[ in ] aRequestedPdofValues    Reference to the matrix of requested pdof values
                 * @param[ in ] aIsMaster             enum for master or slave
                 */
                void get_my_pdof_values(
                        const moris::Cell< Matrix< DDRMat > >  & aPdofValues,
                        const moris::Cell< enum Dof_Type >     & aRequestedDofTypes,
                        Cell< Cell< Matrix< DDRMat > > > & aRequestedPdofValues,
                        const mtk::Master_Slave                  aIsMaster = mtk::Master_Slave::MASTER );

                //-------------------------------------------------------------------------------------------------

                void reshape_pdof_values(
                        const Cell< Matrix< DDRMat > > & aPdofValues,
                        Matrix< DDRMat >         & aReshapedPdofValues );

                //-------------------------------------------------------------------------------------------------

                void set_vector_entry_number_of_pdof();

                //-------------------------------------------------------------------------------------------------

                void get_egn_obj_jacobian( Matrix< DDRMat > & aEqnObjMatrix );

                //-------------------------------------------------------------------------------------------------

                void get_equation_obj_residual( Cell< Matrix< DDRMat > > & aEqnObjRHS );

                //-------------------------------------------------------------------------------------------------

                void add_staggered_contribution_to_residual( Matrix< DDRMat > & aElementResidual );

                //-------------------------------------------------------------------------------------------------
                //            void get_equation_obj_dof_ids( Matrix< DDSMat > & aEqnObjAdofId )
                //            {
                //                aEqnObjAdofId = mUniqueAdofList;
                //            };

                //-------------------------------------------------------------------------------------------------
                void get_equation_obj_dof_ids( Matrix< DDSMat > & aEqnObjAdofId );

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

                //------------------------------------------------------------------------------
                /**
                 * compute dRdp
                 */
                virtual void compute_dRdp()
                {
                    MORIS_ERROR( false, "Equation_Object::compute_dRdp - not implemented in msi." );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp with finite difference
                 */
                virtual void compute_dQIdp_explicit()
                {
                    MORIS_ERROR( false, "Equation_Object::compute_dQIdp_explicit - not implemented in msi." );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdp
                 */
                virtual void compute_dQIdp()
                {
                    MORIS_ERROR( false, "Equation_Object::compute_dQIdp - not implemented in msi." );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute dQIdu
                 */
                virtual void compute_dQIdu()
                {
                    MORIS_ERROR( false, "Equation_Object::compute_dQIdu - not implemented in msi." );
                };

                //------------------------------------------------------------------------------
                /**
                 * compute QI
                 */
                virtual void compute_QI()
                {
                    MORIS_ERROR( false, "Equation_Object::compute_QI - not implemented in msi." );
                };

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

                ////-------------------------------------------------------------------------------------------------
                //            /**
                //             * set the list of side ordinals
                //             */
                //            void set_list_of_side_ordinals( const Matrix< IndexMat > & aListOfSideOrdinals )
                //            {
                //                mListOfSideOrdinals = aListOfSideOrdinals;
                //            }
                //
                ////-------------------------------------------------------------------------------------------------
                //            /**
                //             * set the list of time ordinals
                //             */
                //            void set_list_of_time_ordinals( const Matrix< IndexMat > & aListOfTimeOrdinals )
                //            {
                //                mListOfTimeOrdinals = aListOfTimeOrdinals;
                //            }

                //-------------------------------------------------------------------------------------------------

                /**
                 * how many nodes are connected to this element
                 */
                uint get_num_nodes() const
                {
                    return mNodeObj( 0 ).size();
                }

                //-------------------------------------------------------------------------------------------------

                virtual moris::real get_element_nodal_pdof_value( moris_index   aVertexIndex,
                        moris::Cell< MSI::Dof_Type > aDofType )
                {
                    MORIS_ERROR( false, "Equation_Object::get_element_nodal_pdof_value - this function does nothing");
                    return 0.0;
                }

                //-------------------------------------------------------------------------------------------------

                virtual void set_visualization_cluster( const mtk::Cluster * aVisMeshCluster )
                {
                    MORIS_ASSERT( false, "set_visualization_cluster(), not implemented for base class." );
                }

                //-------------------------------------------------------------------------------------------------
                virtual void compute_quantity_of_interest( const uint            aMeshIndex,
                        enum vis::Output_Type aOutputType,
                        enum vis::Field_Type    aFieldType)
                {
                    MORIS_ASSERT( false, "compute_quantity_of_interest(), not implemented for base class." );
                }

                //-------------------------------------------------------------------------------------------------
        };
    }
}

#endif /* SRC_FEM_CL_EQUATION_OBJECT_HPP_ */
