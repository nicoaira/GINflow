cmake_minimum_required(VERSION 3.24)

option(PQ_CAGRA_CPU_ONLY "Build CPU search/load only (no nvcc)" OFF)

if(PQ_CAGRA_CPU_ONLY)
  project(pq_cagra_adc LANGUAGES CXX)
else()
  project(pq_cagra_adc LANGUAGES CXX CUDA)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT PQ_CAGRA_CPU_ONLY)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 70-real 75-real 80-real 86-real)
  endif()
endif()

if(NOT DEFINED Python_EXECUTABLE AND DEFINED ENV{Python_EXECUTABLE})
  set(Python_EXECUTABLE $ENV{Python_EXECUTABLE})
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(OpenMP)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
  OUTPUT_VARIABLE pybind11_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY
)
find_package(pybind11 CONFIG REQUIRED)

if(PQ_CAGRA_CPU_ONLY)
  pybind11_add_module(pq_cagra_adc_ext MODULE
    src/index_cpu.cpp
    src/index_io.cpp
    src/search_cpu.cpp
    src/bindings.cpp
  )
else()
  find_package(CUDAToolkit REQUIRED)
  pybind11_add_module(pq_cagra_adc_ext MODULE
    src/cagra_pq_adc.cu
    src/index_io.cpp
    src/search_cpu.cpp
    src/bindings.cpp
  )
  target_link_libraries(pq_cagra_adc_ext PRIVATE CUDA::cudart_static)
  set_target_properties(pq_cagra_adc_ext PROPERTIES
    CUDA_SEPARABLE_COMPILATION OFF
  )
  target_compile_options(pq_cagra_adc_ext PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
  )
endif()

target_include_directories(pq_cagra_adc_ext PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
if(OpenMP_CXX_FOUND)
  target_link_libraries(pq_cagra_adc_ext PRIVATE OpenMP::OpenMP_CXX)
endif()

execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('purelib'))"
  OUTPUT_VARIABLE PY_PURELIB
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY
)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(PQ_CAGRA_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/pq_cagra_adc")
else()
  set(PQ_CAGRA_INSTALL_DIR "${PY_PURELIB}/pq_cagra_adc")
endif()

add_custom_command(
  TARGET pq_cagra_adc_ext
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/python/pq_cagra_adc"
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:pq_cagra_adc_ext> "${CMAKE_CURRENT_SOURCE_DIR}/python/pq_cagra_adc/"
)

install(TARGETS pq_cagra_adc_ext LIBRARY DESTINATION "${PQ_CAGRA_INSTALL_DIR}")
install(FILES python/pq_cagra_adc/__init__.py DESTINATION "${PQ_CAGRA_INSTALL_DIR}")
