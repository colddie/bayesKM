project(libtpcimgio)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpcimgio SHARED ${TPC_USE_SOURCE} )

target_include_directories(libtpcimgio PRIVATE ../include)