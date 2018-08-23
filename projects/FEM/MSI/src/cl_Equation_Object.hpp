/*
 * cl_Equation_Object.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_EQUATION_OBJECT_HPP_
#define SRC_FEM_CL_EQUATION_OBJECT_HPP_

#include <memory>

#include "linalg.hpp"
#include "cl_Pdof_Host.hpp"
//#include "cl_FEM_Element.hpp"
//#include "cl_FEM_IWG.hpp"
#include "cl_MSI_Node.hpp"

namespace moris
{
//FIXME will be removed soon
class Linear_Solver;
    namespace MSI
    {
    class Pdof_Host;
    class Equation_Object
    {

    protected:
        moris::Cell< msi::Node * >   mNodeObj;
        moris::uint                   mElementID;
        moris::Cell< Pdof_Host * >    mMyPdofHosts;             // Pointer to the pdof hosts of this equation object

        moris::Cell< enum Dof_Type > mEqnObjDofTypeList;
        moris::Mat< moris::uint >    mTimeSteps;

        moris::Cell< Pdof* >         mFreePdofs;

        moris::Mat< moris::sint >               mUniqueAdofList; // Unique adof list for this equation object
        moris::map < moris::uint, moris::uint > mUniqueAdofMap;  // FIXME replace this map with an MAT. is basically used like a map right now

        moris::Mat< moris::real > mResidual;
        moris::Mat< moris::real > mJacobian;

        moris::Mat< moris::real > mPdofValues;

        //moris::fem::Element* mElement = nullptr;

        // Integrationorder for dof types

        // dof types eg Temp

        //FIXME will be deleted soon. just for testing
        std::shared_ptr< Linear_Solver > mLin;

//-------------------------------------------------------------------------------------------------
    public:
//-------------------------------------------------------------------------------------------------
        Equation_Object()
        {
//            mDofType1.resize( 2, Dof_Type::TEMP );
//            mDofType1( 1 ) = Dof_Type::UX;
        };

//-------------------------------------------------------------------------------------------------
        Equation_Object( const moris::Cell< msi::Node* > & aNodeObjs ) : mNodeObj( aNodeObjs )
        {
            mTimeSteps.resize( 1, 1 );
            mTimeSteps( 0, 0 ) = 0;
        };

    //-------------------------------------------------------------------------------------------------

        ~Equation_Object()
        {
            /*// delete element pointer if it was created
            if ( mElement != NULL )
            {
                delete mElement;
            } */
        };


    //-------------------------------------------------------------------------------------------------
        void get_dof_types( moris::Cell< enum Dof_Type > &  aDofType )
        {
            aDofType = mEqnObjDofTypeList;
        }

    //-------------------------------------------------------------------------------------------------
        const moris::uint get_num_pdof_hosts()
        {
            // FIXME flag used HMR nodes

            // Number of potential pdof hosts based on the number of nodes // Fixme add elements and ghosts
            moris::uint tMyPdofHosts = mNodeObj.size();
            return tMyPdofHosts;
        }

    //-------------------------------------------------------------------------------------------------
        const moris::uint get_max_pdof_hosts_ind()
        {
            moris::luint tMaxPdofHostsInd = 0;

            for ( moris::uint Ii=0; Ii < mNodeObj.size(); Ii++ )
            {
                tMaxPdofHostsInd = std::max( tMaxPdofHostsInd, mNodeObj( Ii )->get_id() );
            }
            return ( moris::uint ) tMaxPdofHostsInd;
        }


    //-------------------------------------------------------------------------------------------------
        void create_my_pdof_hosts( const moris::uint                    aNumUsedDofTypes,
                                   const moris::Mat< moris::sint >    & aPdofTypeMap,
                                         moris::Cell< Pdof_Host * >   & aPdofHostList)
        {
            // Determine size of list containing this equations objects pdof hosts
            moris::uint tNumMyPdofHosts = mNodeObj.size();        //Fixme Add ghost and element numbers

            // Resize list containing this equations objects pdof hosts
            mMyPdofHosts.resize( tNumMyPdofHosts, nullptr );

            // Loop over all nodes of this element, creating new pdof hosts if not existing yet.
            for ( moris::uint Ii=0; Ii < mNodeObj.size(); Ii++ )
            {
                // Save node id Ii in temporary variable for more clarity.
                auto tNodeID = mNodeObj( Ii )->get_id();

                // check if pdof host corresponding to this node exists.
                if ( aPdofHostList( tNodeID ) == NULL)
                {
                    // If node does not exist, create new pdof host.
                    aPdofHostList( tNodeID ) = new Pdof_Host( aNumUsedDofTypes, mNodeObj( Ii ) );
                }
                else
                {
                    // Fixme add else
                }

                // Add pointer to pdof host to the list containing this equation objects pdof hosts.
                mMyPdofHosts( Ii ) = aPdofHostList( tNodeID );

                // FIXME rewrite this function
                for ( moris::uint Ik=0; Ik < mEqnObjDofTypeList.size(); Ik++ )
                {
                    mMyPdofHosts( Ii )->set_pdof_type( mEqnObjDofTypeList( Ik ), mTimeSteps, aNumUsedDofTypes, aPdofTypeMap );
                }
            }
            // Fixme add element

           // FIXME return pointer to pdofs
        }

    //-------------------------------------------------------------------------------------------------
        void create_my_pdof_list()
        {
            // Get number of pdof hosts corresponding to this equation object
            moris::uint tNumMyPdofHosts = mMyPdofHosts.size();

            // Loop over all pdof hosts and get their number of (free) pdofs
            moris::uint tNumMyFreePdofs = 0;
            for ( moris::uint Ii=0; Ii < tNumMyPdofHosts; Ii++ )
            {
                tNumMyFreePdofs = tNumMyFreePdofs + mMyPdofHosts( Ii )->get_num_pdofs();
            }

            // Set size of vector containing this equation objects free pdofs.
            mFreePdofs.resize( tNumMyFreePdofs );

            // FIXME use mempointer instead of loops (?)
            // Loop over all pdof hosts getting pointers to their pdofs. Put these pointers into this equation objects pdof list
            moris::uint tCounter = 0;
            for ( moris::uint Ik=0; Ik < tNumMyPdofHosts; Ik++ )
            {
                // Loop over all pdof types
                for ( moris::uint Ij=0; Ij < ( mMyPdofHosts( Ik )->get_pdof_hosts_pdof_list() ).size(); Ij++ )
                {
                    // Loop over all pdof types times // FIXME
                    for ( moris::uint Ia=0; Ia < ( mMyPdofHosts( Ik )->get_pdof_hosts_pdof_list() )( Ij ).size(); Ia++ )
                    {
                        // add pdof pointer to this equation objects pdof pointer list
                        mFreePdofs( tCounter ) = ( mMyPdofHosts( Ik )->get_pdof_hosts_pdof_list() )( Ij )( Ia );
                        tCounter = tCounter + 1;
                    }
                }
            }
        };

    //-------------------------------------------------------------------------------------------------
        void create_my_list_of_adof_ids()
        {
            // Get MAX number of pdofs for this equation object
            moris::uint tNumMyPdofs = mFreePdofs.size();

            // Loop over all pdofs to count their adofs
            moris::uint tNumMyAdofs = 0;
            for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
            {
                // Get Number of adofs cooresponding to this pdof
                moris::uint tNumAdofForThisPdof = ( mFreePdofs( Ij )->mAdofIds ).length();
                tNumMyAdofs = tNumMyAdofs + tNumAdofForThisPdof;
            }

            // Temporary matrix for adofs Ids
            moris::Mat< moris::sint > tNonUniqueAdofIds( tNumMyAdofs, 1 );

            moris::uint tAdofPosCounter = 0;

            // Loop over all pdofs to get their adofs and put them into a unique list
            for ( moris::uint Ij=0; Ij < tNumMyPdofs; Ij++ )
            {
                tNonUniqueAdofIds ( {tAdofPosCounter, tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).length() -1 }, { 0, 0} ) = mFreePdofs( Ij )->mAdofIds.data();

                // Add number if these adofs to number of assembled adofs
                tAdofPosCounter =tAdofPosCounter + ( mFreePdofs( Ij )->mAdofIds ).length();
            }
            // make list unique
            mUniqueAdofList = moris::unique( tNonUniqueAdofIds );
        };

    //-------------------------------------------------------------------------------------------------
        void set_unique_adof_map()
        {
            //Get number of unique adofs of this equation object
            moris::uint tNumUniqueAdofs = mUniqueAdofList.length();

            // Loop over all unique adofs of this equation object
            for ( moris::uint Ii = 0; Ii < tNumUniqueAdofs; Ii++ )
            {
                mUniqueAdofMap[ mUniqueAdofList( Ii, 0 ) ] = Ii;
            }
        };

    //-------------------------------------------------------------------------------------------------
        // FIXME return map not as a copy. perhaps as input
        void build_PADofMap( moris::Mat< moris::real > & aPADofMap )
        {
            //Get number of unique adofs of this equation object
            moris::uint tNumUniqueAdofs = mUniqueAdofList.length();

            // Get MAX number of pdofs for this equation object
            moris::uint tNumMyPdofs = mFreePdofs.size();

            aPADofMap.set_size( tNumMyPdofs, tNumUniqueAdofs, 0.0 );

            // Loop over all pdofs of this equation object
            for ( moris::uint Ii = 0; Ii < tNumMyPdofs; Ii++ )
            {
                // Loop over all adof Ids of this pdof
                for ( moris::uint Ik = 0; Ik < mFreePdofs( Ii )->mAdofIds.length(); Ik++ )
                {
                    // Getting tPADofMap column entry for the corrsponding value
                    moris::uint tColumnPos = mUniqueAdofMap[ mFreePdofs( Ii )->mAdofIds( Ik, 0 ) ];

                    // Insert value into pdof-adof-map
                    aPADofMap( Ii, tColumnPos ) = ( mFreePdofs( Ii )->mTmatrix)( Ik, 0 );
                }
            }
        };

        //-------------------------------------------------------------------------------------------------
        void get_egn_obj_jacobian( moris::Mat< moris::real > & aEqnObjMatrix )
        {
            moris::Mat< moris::real> tTMatrix;

            this->build_PADofMap( tTMatrix );

            aEqnObjMatrix = trans( tTMatrix )* mJacobian *  tTMatrix ;
        };

        //-------------------------------------------------------------------------------------------------
        void get_equation_obj_residual( moris::Mat< moris::real > & aEqnObjRHS )
        {
            moris::Mat< moris::real> tTMatrix;

            this->build_PADofMap( tTMatrix );

            aEqnObjRHS = trans( tTMatrix ) * mResidual;
        };

        void get_equation_obj_dof_ids( moris::Mat< int > & aEqnObjAdofId )
        {
            aEqnObjAdofId = mUniqueAdofList;

        };

        void get_pdof_values( Mat < real > & aValues );

        //FIXME will be deleted soon
        void set_solver( std::shared_ptr< Linear_Solver > aLin);
    };
    }
}

#endif /* SRC_FEM_CL_EQUATION_OBJECT_HPP_ */
