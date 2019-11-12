project(libtpcimgp)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpcimgp SHARED ${TPC_USE_SOURCE})

target_include_directories(libtpcimgp PRIVATE ../include)