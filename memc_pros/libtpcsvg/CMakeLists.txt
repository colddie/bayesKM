project(libtpcsvg)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpcsvg SHARED ${TPC_USE_SOURCE})

target_include_directories(libtpcsvg PRIVATE ../include)