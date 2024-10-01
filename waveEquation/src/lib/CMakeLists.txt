include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${CMAKE_CURRENT_SOURCE_DIR}/../extern/sepVector/lib)
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/../extern/sepVector/external/hypercube/lib)

# LIST THE SOURCES AND INCLUDES FOR OUR  AWESOME LIBRARY
set(SOURCES wave.cc)
set(INCLUDES wave.h)

add_custom_command(
  OUTPUT kernel.o
  COMMAND ispc --target=avx2-i32x8 ${CMAKE_CURRENT_SOURCE_DIR}/kernel.ispc -o
          kernel.o
  DEPENDS kernel.ispc)

add_library(projWave kernel.o ${SOURCES} ${INCLUDES})
install(TARGETS projWave DESTINATION lib)
