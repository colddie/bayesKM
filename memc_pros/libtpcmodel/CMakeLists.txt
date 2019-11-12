project(libtpcmodel)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpcmodel SHARED ${TPC_USE_SOURCE})

target_include_directories(libtpcmodel PRIVATE ../include)