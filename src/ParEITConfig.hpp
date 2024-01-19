#ifndef __PAREITCONFIG_HPP__
#define __PAREITCONFIG_HPP__

#define EIT_EPSILON 1E-5
#define EIT_EPSILON_DIC 1E-8
#define EIT_ITMAX 1E4
#define EIT_DIM 2
#define EIT_ELEC 4
#define EPSILON 1.0E-12
#define EIT_NAN std::numeric_limits<real_t>::quiet_NaN()
#define EIT_INVALID_INT std::numeric_limits<int_t>::min()
#define EIT_DEBUG 1
#define EIT_PI 3.14159265358979323846 /* pi */

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#define EIT_MPI_REAL MPI_DOUBLE
using real_t = double;
#define EIT_MPI_INT MPI_INT
using int_t = int;
#define EIT_MPI_UINT MPI_UNSIGNED
using uint_t = unsigned int;
using uint   = uint_t;

#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using ElectrodeThetas = std::array<real_t, 2>;

template <class _Type>
struct EitDataType
{
    static MPI_Datatype
    get()
    {
        return MPI_DATATYPE_NULL;
    }
};

#define EIT_DECL_DATAT_TYPE( x, y )                                                                \
    template <>                                                                                    \
    struct EitDataType<x>                                                                          \
    {                                                                                              \
        static MPI_Datatype                                                                        \
        get()                                                                                      \
        {                                                                                          \
            return y;                                                                              \
        }                                                                                          \
    }

EIT_DECL_DATAT_TYPE( real_t, EIT_MPI_REAL );
EIT_DECL_DATAT_TYPE( uint_t, EIT_MPI_UINT );
EIT_DECL_DATAT_TYPE( int_t, EIT_MPI_INT );

namespace parallel
{
#define EIT_MASTER_PROC 0
#define EIT_PROC_RANK parallel::rank
#define EIT_PROC_SIZE parallel::size

extern int rank;
extern int size;
} // namespace parallel

#define EIT_DECL_PARALLEL_DEFINITIONS                                                              \
    namespace parallel                                                                             \
    {                                                                                              \
    int rank;                                                                                      \
    int size;                                                                                      \
    }


#ifndef NDEBUG
    #define ERROR(code) code
#else
    void setErr();
    bool isErr();
    #define ERROR(code)     \
        setErr();           \
        return;
#endif

#endif /* __PAREITCONFIG_HPP__ */