/*
 * cl_Dof_Manager.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_DOF_MANAGER_HPP_
#define SRC_FEM_CL_DOF_MANAGER_HPP_

#include "cl_Cell.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_Equation_Object.hpp"

#include "fn_sum.hpp"

namespace moris
{
namespace MSI
{
class Pdof_Host;
class Dof_Manager
{
private:
    moris::Cell < Pdof_Host * >  mPdofHostList;           // List of all pdof hosts
    moris::Cell < Adof * >       mAdofList;               // List of all adofs
    moris::Cell < Adof * >       mAdofListOwned;          // List of all owned adofs
    //moris::Cell < Adof * >       mAdofListShared;

    moris::Cell< enum Dof_Type > mPdofTypeList;            // List containing all used unique dof types.
    moris::Mat< moris::sint >    mPdofTypeMap;             // Map which maps the unique dof types onto consecutive values.
    moris::Mat< moris::uint >    mPdofHostTimeLevelList;   // List containing the number of time levels per dof type.
    moris::Mat< moris::uint >    mCommTable;               // Communication table. As and input from the model.

    //-----------------------------------------------------------------------------------------------------------
    const moris::uint initialize_max_number_of_possible_pdof_hosts( moris::Cell < Equation_Object* > & aListEqnObj );

    void communicate_dof_types( moris::Cell< enum Dof_Type > & aPdofTypeList );

    void create_dof_type_map();

    void initialize_pdof_host_time_level_list();

    void communicate_time_list( moris::Mat< moris::uint > & aTimeLevelList );

    void communicate_check_if_owned_adof_exists( moris::Cell< moris::Cell < Adof * > > & tAdofListofTypes );

    const moris::uint communicate_adof_offsets( const moris::uint & aNumOwnedAdofs );

    void communicate_shared_adof_ids(const moris::Cell< moris::Cell < Adof * > > & aAdofListofTypes,
                                           moris::Mat< moris::uint >             & aListSharedAdofIds,
                                           moris::Mat< moris::uint >             & aListSharedAdofPos);


public:
    Dof_Manager()
    {};

    Dof_Manager(       moris::Cell < Equation_Object* > & aListEqnObj,
                 const moris::Mat< moris::uint >          aCommTable ) : mCommTable( aCommTable )
    {
        this->initialize_pdof_type_list( aListEqnObj );

        this->initialize_pdof_host_list( aListEqnObj );

        this->create_adofs();

        this->set_pdof_t_matrix();

        for ( moris::uint Ii=0; Ii < aListEqnObj.size(); Ii++ )
        {
            aListEqnObj( Ii )->create_my_pdof_list();
            aListEqnObj( Ii )->create_my_list_of_adof_ids();

            aListEqnObj( Ii )->set_unique_adof_map();
        }
    };

    ~Dof_Manager();

    void initialize_pdof_type_list( moris::Cell < Equation_Object* > & aListEqnObj );

    void initialize_pdof_host_list( moris::Cell < Equation_Object* > & aListEqnObj );

    void create_adofs();

    void set_pdof_t_matrix();

    const moris::uint get_num_adofs() { return mAdofListOwned.size(); };

    const moris::Mat< int > get_local_adof_ids();

    const moris::Mat< int > get_local_overlapping_adof_ids();

};
}
}

#endif /* SRC_FEM_CL_DOF_MANAGER_HPP_ */
