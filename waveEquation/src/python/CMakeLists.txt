project(hypercube-python)
cmake_minimum_required(VERSION 3.0)

set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

#require python3
#find_package(PythonInterp 3 REQUIRED)

#need the lib directory for include
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../lib")

#need pybind11
#find_package(pybind11 REQUIRED)


#need python include directory
include_directories(${PYTHON_INCLUDE_DIR})

#create pybind11 application
pybind11_add_module(pyAcousticDensity MODULE pyacoustic.cc)
target_link_libraries(pyAcousticDensity PRIVATE  PRIVATE  projWave   sepVector hypercube ${OpenMP_EXE_LINKER_FLAGS} ${TBB_LIBRARIES})

#install the libraries
install( TARGETS pyAcousticDensity DESTINATION lib/python)

#We will wrap the python/C++ interface with a strict python
#install(FILES ProjSmooth.py DESTINATION lib/python)
