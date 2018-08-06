/*
 * cl_Dof_Manager.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_DOF_MANAGER_HPP_
#define SRC_FEM_CL_DOF_MANAGER_HPP_

#include "cl_Cell.hpp"
#include "cl_Equation_Object.hpp"

namespace moris
{
namespace MSI
{
class Pdof_Host;
class Dof_Manager
{
private:
    moris::Cell < Pdof_Host * > mPdofHostList;
    moris::Cell < Adof * > mAdofList;

    moris::uint mMaxNumPdofHosts;

    // List containing all dof types. Perhaps make temporary in future.
    moris::Cell< enum Dof_Type > mPdofTypeList;

public:
    Dof_Manager()
    {
        //this->initialize_max_number_of_possible_pdof_hosts( aListEqnObj );
    };

    Dof_Manager( const moris::uint aNumEquationObj,
                       moris::Cell < Equation_Object* > & aListEqnObj )
    {
        this->initialize_max_number_of_possible_pdof_hosts( aListEqnObj );

        this->initialize_pdof_host_list( aListEqnObj );

        this->create_adofs();

        this->set_pdof_t_matrix();

        for ( moris::uint Ii=0; Ii < aListEqnObj.size(); Ii++ )
        {
            aListEqnObj( Ii )->create_my_list_of_adof_ids();

            aListEqnObj( Ii )->set_unique_adof_map();
        }
    };

    ~Dof_Manager()
    {};

//-----------------------------------------------------------------------------------------------------------
    void initialize_max_number_of_possible_pdof_hosts( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        // Ask how many equation objects
        //std::cout<<"num equation objects: "<<aListEqnObj.size()<<std::endl;
        moris::uint tNumEqnObj = aListEqnObj.size();
        mMaxNumPdofHosts = 0;

        //loop over all equation objects, asking for their number of pdof hosts
        for ( moris::uint Ii=0; Ii < tNumEqnObj; Ii++ )
        {
            mMaxNumPdofHosts += aListEqnObj( Ii )->get_num_pdof_hosts();
        }
    };

//-----------------------------------------------------------------------------------------------------------
    void initialize_pdof_host_list( moris::Cell < Equation_Object* > & aListEqnObj )
    {
        // reserve the maximal needed amount of memory for  a temporary PdofHostList
        moris::Cell < Pdof_Host* > tPdofHostList;

        // FIXME resize list of pdof types
        //enum Dof_Type tDofType = Dof_Type::INITIALIZE_DOF_TYPE;
        mPdofTypeList.resize( 256, Dof_Type::INITIALIZE_DOF_TYPE );

        tPdofHostList.resize( mMaxNumPdofHosts, nullptr );

        // Loop over all equation objects. Ask them to create their pdof hosts. Pdof hosts are stored in tPdofHostList
        for ( moris::uint Ii=0; Ii < aListEqnObj.size(); Ii++ )
        {
            aListEqnObj( Ii )->create_my_pdof_hosts( tPdofHostList, mPdofTypeList );
        }

        // Determine number of Pdof Hosts
        moris::uint tNumPdofHosts = 0;
        for ( moris::uint Ik = 0; Ik < tPdofHostList.size(); Ik++ )
        {
            // If pointer in temporary pdof host list exists add one to number of pdof hosts
            if ( tPdofHostList( Ik ) != NULL)
            {
                 tNumPdofHosts = tNumPdofHosts + 1;
            }
        }

        // Set size of List containing all pdof hosts
        mPdofHostList.resize( tNumPdofHosts );

        moris::uint counter = 0;
        // add pointers to pdof hosts into list of pdof hosts
        for ( moris::uint Ij = 0; Ij < tPdofHostList.size(); Ij++ )
        {
            // If pointer in temporary pdof host list exists add one to number of pdof hosts
            if ( tPdofHostList( Ij ) != NULL)
            {
                mPdofHostList( counter ) = tPdofHostList( Ij );
                counter = counter + 1;
            }
        }
    };

//-----------------------------------------------------------------------------------------------------------
    void create_adofs()
    {
        // Get number of pdoftypes
        moris::uint tNumPdofTypes = 0;
        for ( moris::uint Ik = 0; Ik < 258; Ik++ )
        {
            if ( mPdofTypeList( Ik ) == Dof_Type::INITIALIZE_DOF_TYPE )
            {
                tNumPdofTypes = Ik;
                break;
            }
        }

        // Get number of pdof hosts in pdof host list
        moris::uint tNumPdofHosts = mPdofHostList.size();

        // Get max entry of node adof
        moris::sint tMaxNodeAdofId = 0;
        for ( moris::uint Ik = 0; Ik < tNumPdofHosts; Ik++ )
        {
            tMaxNodeAdofId = std::max( tMaxNodeAdofId, mPdofHostList( Ik )->get_node_obj_ptr()->get_adofs().max() );
        }
        // Add one because c++ is 0 based. ==> List size has to be tMaxNodeAdofId + 1
        tMaxNodeAdofId = tMaxNodeAdofId +1;

        // Create temporary moris::Cell containing lists of temporary adofs
        moris::Cell<moris::Cell < Adof * > > tAdofListofTypes;
        tAdofListofTypes.resize( tNumPdofTypes );

        for ( moris::uint Ik = 0; Ik < tNumPdofTypes; Ik++ )
        {
            tAdofListofTypes( Ik ).resize( tMaxNodeAdofId, nullptr );
        }
        //------------------------------------------------
        // Set temporary adof vector. Size equals size of num pdofs
        //moris::Cell < Adof * > tAdofList( tNumPdofs, nullptr );

        // Loop over all pdof hosts and get the adofs
        for ( moris::uint Ii = 0; Ii < tNumPdofHosts; Ii++ )
        {
            mPdofHostList( Ii )->get_adofs( tAdofListofTypes );
        }

        // Determine number of adofs
        moris::uint tNumAdofs = 0;
        for ( moris::uint Ik = 0; Ik < tAdofListofTypes.size(); Ik++ )
        {
            for ( moris::uint Ia = 0; Ia < tAdofListofTypes( Ik ).size(); Ia++ )
            {
                // If pointer in temporary adof list exists. Add one to number of adofs
                if ( tAdofListofTypes( Ik )( Ia ) != NULL)
                {
                    tNumAdofs = tNumAdofs + 1;
                }
            }
        }

        // Set size of List containing all adofs
        mAdofList.resize( tNumAdofs );

        moris::uint tCounter = 0;
        // add pointers to adofs into list of adofs
        for ( moris::uint Ij = 0; Ij < tAdofListofTypes.size(); Ij++ )
        {
            for ( moris::uint Ib = 0; Ib < tAdofListofTypes( Ij ).size(); Ib++ )
            {
                // If pointer in temporary adofs list exists, add one to number of adofs
                if ( tAdofListofTypes( Ij )( Ib ) != NULL)
                {
                    mAdofList( tCounter ) = tAdofListofTypes( Ij )( Ib ) ;

                    // Set adof Id. Fixme introduce offset for parallel
                    mAdofList( tCounter )->set_adof_id( tCounter );

                    tCounter = tCounter + 1;
                }
            }
        }

        // Tell pdofs to get adof Ids
        for ( moris::uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
        {
            // all pdofs of this pdof host will ask for their adof Ids
            mPdofHostList(Ij)->get_adofs_ids();

            // create unique adof list for this pdof host
            mPdofHostList(Ij)->create_unique_adof_list();

            //mPdofHostList(Ij)->set_unique_adof_map();

        }
    };

    //-----------------------------------------------------------------------------------------------------------
    void set_pdof_t_matrix()
    {
        // Get number of pdof hosts in pdof host list
        moris::uint tNumPdofHosts = mPdofHostList.size();

        // Tell pdofs to get t their t matrix
        for ( moris::uint Ij = 0; Ij < tNumPdofHosts; Ij++ )
        {
            // all pdofs of this pdof host will ask for their t matrix
            mPdofHostList(Ij)->set_t_matrix();
        }
    };

    //-----------------------------------------------------------------------------------------------------------

    const moris::uint get_num_adofs()
    {
        moris::uint tNumAdofs = mAdofList.size();
        return tNumAdofs;
    };

    const moris::Mat< int > get_local_adof_ids()
    {
        Mat< int > tLocalAdofIds ( mAdofList.size(), 1 );
        for ( moris::uint Ij = 0; Ij < mAdofList.size(); Ij++ )
        {
            tLocalAdofIds( Ij, 0 ) = mAdofList( Ij )->get_adof_id();
        }
        return tLocalAdofIds;
    };
};
}
}

#endif /* SRC_FEM_CL_DOF_MANAGER_HPP_ */
