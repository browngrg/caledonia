#ifndef MPITypeTraits_H_
#define MPITypeTraits_H_


#include<mpi.h>

/** \file MPITypeTraits.hpp
 *  \author Gregory Brown
 *  \brief A convenient way to handle templated overloading of MPI functions
 */


// namespace psimag {


/** \class MPITypeTraits
 *  \author Gregory Brown
 *  \brief Traits class that converts c++ type to MPI_TYPE.
 *
 *  This is useful for creating templated instances of massively
 *  parallel calls. Inside the arguments to the MPI calls, you
 *  convert the type with 
 *  \verbatim typename psimag::MPITypeTraits<template_type>::mpitype \endverbatim
 */
template<class T> struct MPITypeTraits          {  };

template<> struct MPITypeTraits<char> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<char>::mpitype=MPI_CHAR; 

template<> struct MPITypeTraits<double> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<double>::mpitype=MPI_DOUBLE; 

template<> struct MPITypeTraits<float> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<float>::mpitype=MPI_FLOAT;  

template<> struct MPITypeTraits<int> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<int>::mpitype=MPI_INT;  

template<> struct MPITypeTraits<long> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<long>::mpitype=MPI_LONG;  

template<> struct MPITypeTraits<long double> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<long double>::mpitype=MPI_LONG_DOUBLE;  

template<> struct MPITypeTraits<short> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<short>::mpitype=MPI_SHORT;  

template<> struct MPITypeTraits<unsigned char> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<unsigned char>::mpitype=MPI_UNSIGNED_CHAR;  

template<> struct MPITypeTraits<unsigned int> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<unsigned int>::mpitype=MPI_UNSIGNED;  

template<> struct MPITypeTraits<unsigned long> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<unsigned long>::mpitype=MPI_UNSIGNED_LONG;  

template<> struct MPITypeTraits<unsigned short> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<unsigned short>::mpitype=MPI_UNSIGNED_SHORT;  

template<> struct MPITypeTraits<bool> { static const MPI_Datatype mpitype; };
const MPI_Datatype MPITypeTraits<bool>::mpitype=MPI_SHORT;  

// }       // namespace psimag

#endif  // MPITypeTraits_H_
