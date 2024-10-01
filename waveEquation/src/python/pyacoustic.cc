#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "wave.h"
namespace py = pybind11;
using namespace SEP;


PYBIND11_MODULE(pyAcousticDensity, clsOps) {
py::class_<acousticDensity, std::shared_ptr<acousticDensity>>(clsOps, "acousticDensity")
 .def(py::init<std::shared_ptr<float3DReg>,std::shared_ptr<float3DReg>,
   const int,const int, const int, std::shared_ptr<float1DReg>,
   const float,
   const float,
   const int,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> 
   >(), "Initialize ")
  .def("prop", (void (acousticDensity::*)(const int))
	 &acousticDensity::prop, "Propagate nt steps")
   .def("getWavefield",(std::shared_ptr<float3DReg> (acousticDensity::*)())
		   &acousticDensity::getWavefield,"Get stored Wavefield");
py::class_<acousticDensityISPC, std::shared_ptr<acousticDensityISPC>>(clsOps, "acousticDensityISPC")
 .def(py::init<std::shared_ptr<float3DReg>,std::shared_ptr<float3DReg>,
   const int,const int, const int, std::shared_ptr<float1DReg>,
   const float,
   const float,
   const int,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> 
   >(), "Initialize ")
  .def("prop", (void (acousticDensity::*)(const int))
	 &acousticDensity::prop, "Propagate nt steps")
   .def("getWavefield",(std::shared_ptr<float3DReg> (acousticDensity::*)())
		   &acousticDensity::getWavefield,"Get stored Wavefield");

py::class_<acousticDensityISPC_P, std::shared_ptr<acousticDensityISPC_P>>(clsOps, "acousticDensityISPC_P")
 .def(py::init<std::shared_ptr<float3DReg>,std::shared_ptr<float3DReg>,
   const int,const int, const int, std::shared_ptr<float1DReg>,
   const float,
   const float,
   const int,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> ,
   std::vector<std::vector<int>> 
   >(), "Initialize ")
  .def("prop", (void (acousticDensity::*)(const int))
	 &acousticDensity::prop, "Propagate nt steps")
   .def("getWavefield",(std::shared_ptr<float3DReg> (acousticDensity::*)())
		   &acousticDensity::getWavefield,"Get stored Wavefield");
}

