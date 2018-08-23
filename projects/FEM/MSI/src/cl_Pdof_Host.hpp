/*
 * cl_Pdof_Host.hpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_PDOF_HOST_HPP_
#define SRC_FEM_CL_PDOF_HOST_HPP_

#include "cl_Dof_Type_Enums.hpp"
#include "cl_Node_Obj.hpp"
#include "cl_Adof.hpp"

#include "fn_unique.hpp" // LNA/src
#include "cl_Map.hpp" // LNA/src

namespace moris
{
    namespace MSI
    {
    struct Pdof
    {
        moris::uint                  mDofTypeIndex;
        moris::uint                  mTimeStepIndex;
        moris::Mat < moris::sint >   mAdofIds;
        moris::Mat < moris::real >   mTmatrix;

        moris::Cell < Adof* >        mAdofPtrList;              //FIXME delete this list after call to get adof ids or replace it
    };

//-------------------------------------------------------------------------------------------------
    class Pdof_Host
    {
    private:
        moris::Mat< moris::uint >               mPdofTypeExist;         // Vector indicates if dof type exists. FIXME replace by bitset
        moris::Cell< moris::Cell< Pdof* > >     mListOfPdofTimePerType; // List of all pdofs per time per dof type

        moris::Mat< moris::uint >               mUniqueAdofList;        // Unique adof list for this pdof host
        moris::map < moris::uint, moris::uint > mUniqueAdofMap;         // FIXME membe r function tio build this map is never called

        mtk::Vertex * mNodeObj;                                           // FIXME replace base class bei FEM node
        moris::luint  mNodeID;

    protected:

       //FIXME Add interpolation order

    public:
        Pdof_Host()
        {
        };

        Pdof_Host( const moris::uint   aNumUsedDofTypes,
                         mtk::Vertex * aNodeObj ) : mNodeObj( aNodeObj )
        {
            mNodeID = mNodeObj->get_id();

            mPdofTypeExist.set_size( aNumUsedDofTypes, 1, 0 );

            // Set size of list to the number of used nodes
            mListOfPdofTimePerType.resize( aNumUsedDofTypes );
        };

        ~Pdof_Host();

        mtk::Vertex * const get_node_obj_ptr()
        {
            return mNodeObj;
        };

        /**
         * @brief Sets a dof type and time to this pdof host. This function is tested by the test [Pdof_host_set_dof_type]
         *
         * @param[in] aDof_Type          Dof type
         * @param[in] aTimeSteps         Time steps.
         * @param[in] aNumUsedDofTypes   Globally used dof types.
         * @param[in] aPdofTypeMap       Paw which related a dof type enum to the index.
         *
         */
        void set_pdof_type( const enum Dof_Type                  aDof_Type,
                            const moris::Mat< moris::uint >    & aTimeSteps,
                            const moris::uint                    aNumUsedDofTypes,
                            const moris::Mat< moris::sint >    & aPdofTypeMap);

        /**
         * @brief Gets the adofs for all the pdofs in this pdof host. This function is tested by the test [Pdof_Host_Get_Adofs]
         *
         * @param[in] aTimeLevelOffsets  Offsets for this doftype and time      FIXME
         * @param[in] aAdofList          List containing all the adofs.
         *
         */
        void get_adofs( const moris::Mat< moris::uint >            & aTimeLevelOffsets,
                              moris::Cell< moris::Cell< Adof * > > & aAdofList );

        /**
         * @brief Gets the adofs Ids for all the pdofs in this pdof host. This function is tested by the test [Pdof_Host_Get_Adofs]
         *
         * @param[in] aTimeLevelOffsets  Offsets for this doftype and time      FIXME
         * @param[in] aAdofList          List containing all the adofs.
         *
         */
        void get_adofs_ids();

        // FIXME member function not used
        void create_unique_adof_list();

        /**
         * @brief Set the t-matrix values for all the pdofs. This function is tested by the test [Pdof_Host_Get_Adofs]
         *
         */
        void set_t_matrix();

//-------------------------------------------------------------------------------------------------
        // FIXME member function not used
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

        const moris::uint get_num_time_levels_of_type( const moris::uint & aDofTypeInd ) { return mListOfPdofTimePerType( aDofTypeInd ).size(); };

        const moris::uint get_num_pdofs();

        moris::Cell< moris::Cell< Pdof* > > & get_pdof_hosts_pdof_list() { return mListOfPdofTimePerType; }

        void set_pointer_to_Tmatrix()
        {};

//-------------------------------------------------------------------------------------------------
        void set_adof_IDs()
        {};

//-------------------------------------------------------------------------------------------------
    };

    }
}



#endif /* SRC_FEM_CL_PDOF_HOST_HPP_ */
