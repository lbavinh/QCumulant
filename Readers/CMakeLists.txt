set(READERS_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/PicoDst/include
  ${CMAKE_SOURCE_DIR}/PicoDst/src
  ${ROOT_INLCUDE_DIRS}
)

include_directories(${READERS_INCLUDE_DIRECTORIES})

set(READERS_INCLUDE_LIBRARIES
  PicoDst
  ${ROOT_LIBRARIES}
)

set(Readers_h_files
  "${CMAKE_CURRENT_SOURCE_DIR}/IReader.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/PicoDstReader.h"
)

set(Readers_cxx_files
  "${CMAKE_CURRENT_SOURCE_DIR}/IReader.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/PicoDstReader.cxx"
)

set(Readers_LinkDef
	"${CMAKE_CURRENT_SOURCE_DIR}/Readers.LinkDef.h"
)

#---Generate dictionary
ROOT_GENERATE_DICTIONARY(G__Readers
  ${Readers_h_files}
  LINKDEF ${Readers_LinkDef}
)

#---Compile library
add_library(Readers SHARED ${Readers_cxx_files} G__Readers.cxx)
target_link_libraries(Readers ${READERS_INCLUDE_LIBRARIES})