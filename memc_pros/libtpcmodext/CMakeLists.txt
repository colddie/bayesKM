project(libtpcmodext)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpcmodext SHARED ${TPC_USE_SOURCE})

target_include_directories(libtpcmodext PRIVATE ../include)