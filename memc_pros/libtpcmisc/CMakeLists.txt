project(libtpcmisc)

file(GLOB TPC_USE_SOURCE
 "*.h"
 "*.c"
)

add_library(libtpcmisc SHARED ${TPC_USE_SOURCE})

target_include_directories(libtpcmisc PRIVATE ../include)