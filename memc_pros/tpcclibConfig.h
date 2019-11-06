// Enable POSIX extensions
#define _POSIX_C_SOURCE 200809L
// Enable GNU extensions
#define _GNU_SOURCE
// the configured options and settings for tpcclib
#define tpcclib_VERSION_MAJOR 0
#define tpcclib_VERSION_MINOR 7
#define tpcclib_VERSION_PATCH 1
#define tpcclib_COPYRIGHT "(c) 2019 by Turku PET Centre"
// Are we using MinGW compiler on Windows
/* #undef MINGW */
// Do certain functions exist?
#define HAVE_STRCASESTR
#define HAVE_STRDUP
#define HAVE_STRNDUP
/* #undef HAVE_BUILTIN_STRNDUP */
#define HAVE_STRNLEN
/* #undef HAVE_STRLCAT */
/* #undef HAVE_STRLCPY */
#define HAVE_TIMEGM
#define HAVE_GMTIME_R
/* #undef HAVE_GMTIME_S */
#define HAVE_LOCALTIME_R
#define HAVE_STRPTIME
#define HAVE_TIMESPEC_GET
#define HAVE_CLOCK_GETTIME
#define HAVE_GETTIMEOFDAY
#define HAVE_GETPID
#define HAVE_GET_CURRENT_DIR_NAME
/* #undef HAVE__MKDIR */
// Do certain struct members exist?
/* #undef HAVE_TM_GMTOFF */
#define HAVE_TIME_WITH_SYS_TIME
/* #undef HAVE_ST_BIRTHTIME */
// Do we have direct.h?
/* #undef HAVE_DIRECT_H */
// Do we have omp.h (OpenMP)?
#define HAVE_OMP_H

// By default, structs are aligned; uncomment to pack
// TPCCLIB is written so that both should be fine, but implementation of
// packing seems to be variable in different compilers, compiler versions, 
// and platforms, and thus not safe.
//#pragma pack(1)
