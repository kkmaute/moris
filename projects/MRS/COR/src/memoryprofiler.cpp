/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * memoryprofiler.hpp
 *
 */

#include <cstdio>
#include <cstdint>

// define static varibales to store the file handle and filename
static FILE*       fp       = nullptr;
static const char* filename = "memeory_file.log";

// initialize the file handle and filename 
// attribute constructor is used to call this function before main
__attribute__( ( constructor ) ) void
Initialize_memeory_profiler()
{
    //printf( "Initializing.\n" );
    fp = fopen( filename, "w" );
}

// c linkages are used to avoid name mangling
extern "C" {

// declare the real malloc and free functions
void* __real_malloc( size_t size );

//------------------------------------------------------------------------------

// declare the wrapper malloc and free functions , custom malloc implementation
void*
__wrap_malloc( size_t size )
{
    // get the pointer address from the real malloc
    void* ptr = __real_malloc( size );
    
    // if fp has been assigned
    if ( fp != NULL )
    {    // write the malloc call to a file
        fprintf( fp, "M, %lu, %lu\n", (uintptr_t)ptr, size );
    }
    
    //or write directly it to the log file
    //printf( "M, %lu, %lu\n", (uintptr_t)ptr, size );

    // return the pointer address
    return ptr;
}

//------------------------------------------------------------------------------

// declare the real free function
void __real_free( void* ptr );

// define the custom free function
void
__wrap_free( void* ptr )
{
    // Your custom code here.
    if ( fp != NULL )
    { 
        fprintf(fp, "F, %lu\n", (uintptr_t)ptr);
    }

    //or write directly it to the log file
    //printf( "F, %lu\n", (uintptr_t)ptr );

    // Call the real free.
    __real_free( ptr );
}

//------------------------------------------------------------------------------

// declare one of the the memeory allocation function arma uses 
int __real_posix_memalign( void** memptr, size_t alignment, size_t size );

//------------------------------------------------------------------------------

// the custoom version of the posix_memalign function
int __wrap_posix_memalign( void** memptr, size_t alignment, size_t size )
{
    // call the real function and store the return value 
    int tStatus = __real_posix_memalign( memptr, alignment, size );

    if ( fp != NULL )
    { 
        fprintf( fp, "M, %lu, %lu\n", (uintptr_t)(*memptr), size );
    }

    // Your custom code here.
    //printf( "M, %lu, %lu\n", (uintptr_t)(&memptr), size );

    // return the status of the real function
    return tStatus;
}

//------------------------------------------------------------------------------

}

// overload the new and delete operators
void*
operator new( std::size_t size )
{
    return __wrap_malloc( size );
}

//------------------------------------------------------------------------------
void
operator delete( void* ptr ) noexcept
{
    __wrap_free( ptr );
}

//------------------------------------------------------------------------------

void*
operator new[]( std::size_t size )
{
    return __wrap_malloc( size );
}

//------------------------------------------------------------------------------

void
operator delete[]( void* ptr ) noexcept
{
    __wrap_free( ptr );
}
