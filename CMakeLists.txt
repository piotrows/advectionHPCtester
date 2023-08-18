cmake_minimum_required (VERSION 3.14.0)
project(advectionHPCtester Fortran )
enable_language (Fortran)

set(CMAKE_Fortran_FLAGS "-g ")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(optsprecsp -DFLOAT_PRECISION=4) 
set(optsprecdp -DFLOAT_PRECISION=8) 

option(ENABLE_GPU "Configure for NVIDIA GPU" OFF)

if(ENABLE_GPU)
find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
  message(STATUS "Building with GPU support")
  enable_language (CUDA)
  find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
  find_library(LIBPNETCDF libpnetcdf PATHS "/usr/lib/x86_64-linux-gnu")
  set(optsgpu -DCUDACODE -DTIMERSCUDA  )
else()
  message(STATUS "Building without GPU support")
  set(optsgpu -DNOCUDACODE   )
endif()

option(ENABLE_MPI "Configure for MPI parallelization" OFF)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "Building with MPI")
  if(ENABLE_GPU)
    set(optsmpi  -DPUREMPI -DMPI90 )
  else()
    set(optsmpi  -DPUREMPI -DMPI90 -DTIMERSMPI )
  endif()
  message (STATUS ${MPI_Fortran_INCLUDE_PATH})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
else()
  message(STATUS "Building without MPI")
  if(ENABLE_GPU)
  else()
    set(optsmpi -DTIMERSCPU    )
  endif()
endif()


option(ENABLE_PNETCDF "Configure PNETCDF output" OFF)

if(ENABLE_PNETCDF)
 message(STATUS "Building with PNETCDF support")
 set(optspnet  -DPNETCDF -lpnetcdf -I$ENV{PNETCDF_INCLUDE} -L$ENV{PNETCDF_LIB} )
 TARGET_LINK_LIBRARIES(advection common)
else()
 message(STATUS "Building without PNETCDF support")
 set(optspnet   )
endif()

option(ENABLE_TESTING "Compute solution norms" ON)

if(ENABLE_TESTING)
 message(STATUS "Building with TESTING of the solutions ")
 set(optstest   -DTESTING )
else()
        message(STATUS "Building without TESTING support")
 set(optstest "")
endif()

option(ENABLE_FAST "Optimize - ON, Debug - OFF"  ON)

if(ENABLE_FAST)

	#  message(STATUS ${CMAKE_Fortran_COMPILER_ID})
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(optims  -Ofast -cpp -g -fbacktrace )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(optims -fpp  -O3 -xHost -unroll-aggressive -qopt-prefetch -align commons )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(optims  -O2 -Mpreprocess )
  endif()

else()

  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
       set(optims  "-ffree-form -fimplicit-none -cpp"  )
       set(debugs  "-g -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow  -finit-real=snan"  )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
       set(optims  -O1 -fpp)
       set(debugs  -g -traceback -check bounds -ftrapuv -init=snan -implicitnone -fpe0 -ftz )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
       set(optims  -O1 ) #-Mpreprocess )
       set(debugs  -C traceback -Mbounds -Mchkfpstk -Kieee -Mlist -Mpreprocess )
  endif()

endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debugs}")
set(CMAKE_Fortran_FLAGS_TMP  "${optims} ${optsmpi} ${optsgpu} ${optpnet} ${optstest}")
if(ENABLE_GPU)
set(CMAKE_Fortran_FLAGS_TMP "${CMAKE_Fortran_FLAGS_TMP} -cuda")
endif()

foreach(arg ${CMAKE_Fortran_FLAGS_TMP})
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${arg}")
endforeach(arg ${CMAKE_Fortran_FLAGS_TMP})
message(STATUS "Compiling with flags: ${CMAKE_Fortran_FLAGS} ")

set(common_sources ${common_sources}
precisions.f90
bconditions.f90
epsilons.f90
eulag_lib_parameters.f90
noise.f90
argsparser/argsparser.f90)

set(common_sources_preprocess ${common_sources_preprocess}
mpi_parallel.f90
eulag_datafields_modular.f90
diagnostics.f90
diagutils.f90
eulag_moist_datafields.f90
scratch_datafields_modular.f90
geometry_datafields_modular.f90
geometry.f90
filters.f90
)

set(advection_alg_sources_preprocess ${advection_alg_sources_preprocess}
upwind3d-gpubc.f90
antidiff3d-gauge-gpubc.f90
antidiff3d-standard-gpubc.f90
mpdata3d-rho-gauge-gpubc.f90
velprd_driver.f90
module_velprd.f90
)

set(advection_drv_sources_preprocess ${advection_drv_sources_preprocess}
advec.f90
advec_init_finalize.f90
advection_CEtimestep_twostep_GPUBC.f90
)
set(advection_tst_sources_preprocess ${advection_tst_sources_preprocess}
advection_test_setup.f90
)

function(DEFINE_PREPROCESS_SED subdir_arg destdir_arg files_arg custopt_arg)
set(subdir ${subdir_arg})
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/${destdir_arg}/${subdir}") 
set(localdir "${PROJECT_SOURCE_DIR}/${subdir}")
set(destdir  "${destdir_arg}/${subdir}")
foreach( file_i ${files_arg})
        message (STATUS "define preprocessing for file: ${destdir}/${file_i}")
        add_custom_command( 
        COMMAND gcc -E -P ${optsmpi} ${optsgpu} ${optstest} ${custopt_arg} -cpp  ${localdir}/${file_i}
        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${destdir}/${file_i}.pp.f90"
        DEPENDS  ${localdir}/${file_i}
        OUTPUT  ${destdir}/${file_i}.pp.f90
        COMMENT "Preprocessing ${localdir}/${file_i} to ${destdir}/${file_i}.pp.f90"
        VERBATIM
)
endforeach( file_i )
endfunction(DEFINE_PROCESS_SED)

function(GENERATE_BUILD buildtype_arg specopts_arg driverfile_arg) 
set(buildtype ${buildtype_arg})
set(specopts  ${specopts_arg} )
set(driverfile ${driverfile_arg} )
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/common) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/common/argsparser) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/advection) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/include) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/lib) 
set(advection_sources )  

define_preprocess_sed("common" "${buildtype}" "${common_sources}"            "${specopts}")
define_preprocess_sed("common" "${buildtype}" "${common_sources_preprocess}" "${specopts}")
define_preprocess_sed("advection/src_algorithms" "${buildtype}" "${advection_alg_sources_preprocess}" "${specopts}")

define_preprocess_sed("advection/src_drivers"    "${buildtype}" "${advection_drv_sources_preprocess}" "${specopts}")
define_preprocess_sed("advection/src_testdata"   "${buildtype}" "${advection_tst_sources_preprocess}" "${specopts}")

list(TRANSFORM advection_alg_sources_preprocess PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/advection/src_algorithms/)
list(TRANSFORM advection_alg_sources_preprocess APPEND .pp.f90)
list(TRANSFORM advection_drv_sources_preprocess PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/advection/src_drivers/)
list(TRANSFORM advection_drv_sources_preprocess APPEND .pp.f90)
list(TRANSFORM advection_tst_sources_preprocess PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/advection/src_testdata/)
list(TRANSFORM advection_tst_sources_preprocess APPEND .pp.f90)
list(APPEND advection_sources ${advection_alg_sources_preprocess} )
list(APPEND advection_sources ${advection_drv_sources_preprocess} )
list(APPEND advection_sources ${advection_tst_sources_preprocess} )

list(TRANSFORM common_sources            PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/common/)
list(TRANSFORM common_sources            APPEND .pp.f90)
list(TRANSFORM common_sources_preprocess PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/common/)
list(TRANSFORM common_sources_preprocess APPEND .pp.f90)
list(APPEND common_sources ${common_sources_preprocess} )

add_library(common_${buildtype}    ${common_sources})   
add_library(advection_${buildtype} ${advection_sources})
add_dependencies(advection_${buildtype}  common_${buildtype} )
target_link_libraries(advection_${buildtype} PRIVATE ${MPI_Fortran_LIBRARIES} common_${buildtype} )

set(commondir    ${CMAKE_BINARY_DIR}/${buildtype}/include)
set(advectiondir ${CMAKE_BINARY_DIR}/${buildtype}/include)
set_target_properties(   common_${buildtype} PROPERTIES Fortran_MODULE_DIRECTORY ${commondir})
set_target_properties(   common_${buildtype} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/lib)
set_target_properties(advection_${buildtype} PROPERTIES Fortran_MODULE_DIRECTORY ${advectiondir})
set_target_properties(advection_${buildtype} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/lib)
target_include_directories(common_${buildtype} PRIVATE ${commondir} ) 
target_include_directories(advection_${buildtype} PRIVATE ${advectiondir}) 

list(TRANSFORM driverfile PREPEND ${CMAKE_SOURCE_DIR}/advection/src_drivers/)
add_library(advection_interface_${buildtype} SHARED ${driverfile})
add_dependencies(advection_interface_${buildtype}  advection_${buildtype} )
target_include_directories(advection_interface_${buildtype} PRIVATE ${commondir} ) 
target_include_directories(advection_interface_${buildtype} PRIVATE ${advectiondir}) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include) 
set_target_properties(advection_interface_${buildtype} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set_target_properties(advection_interface_${buildtype} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

endfunction(GENERATE_BUILD)

set(buildtype sp)
set(driverfile 
    advec_interface_${buildtype}.f90
    advec_init_interface_${buildtype}.f90
    )
generate_build(${buildtype} ${optsprecsp} ${driverfile}) 
add_link_options( -ldl)
if(ENABLE_GPU)
add_link_options( -cuda -lnvToolsExt)
set(execname_single advectiontest.gpu.${buildtype})
else()
set(execname_single advectiontest.cpu.${buildtype})
endif()
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal_${buildtype}.f90")
add_dependencies(${execname_single} common_${buildtype} advection_${buildtype} )
target_link_libraries(advection_interface_${buildtype} PRIVATE advection_${buildtype})
target_link_libraries(${execname_single}  advection_interface_${buildtype})
target_include_directories(${execname_single} PRIVATE ${CMAKE_BINARY_DIR}/include) 
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)

set(buildtype dp)
set(driverfile 
    advec_interface_${buildtype}.f90
    advec_init_interface_${buildtype}.f90)
generate_build(${buildtype} ${optsprecdp} ${driverfile}) 
if(ENABLE_GPU)
add_link_options( -cuda -lnvToolsExt)
set(execname_single advectiontest.gpu.${buildtype})
else()
set(execname_single advectiontest.cpu.${buildtype})
endif()
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal_${buildtype}.f90")
add_dependencies(${execname_single} common_${buildtype} advection_${buildtype} )
target_link_libraries(advection_interface_${buildtype} PRIVATE  advection_${buildtype})
target_link_libraries(${execname_single}  advection_interface_${buildtype})
target_include_directories(${execname_single} PRIVATE ${CMAKE_BINARY_DIR}/include) # ${solverdir} ${diffusiondir})
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)


set(buildtype uni)
set(buildtypesp sp)
set(buildtypedp dp)
if(ENABLE_GPU)
add_link_options( -cuda -lnvToolsExt)
set(execname_single advectiontest.gpu.${buildtype})
else()
set(execname_single advectiontest.cpu.${buildtype})
endif()
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal.f90")
add_dependencies(${execname_single} advection_interface_${buildtypedp} 
                                    advection_interface_${buildtypesp}  )
target_link_libraries(${execname_single} PRIVATE ${MPI_Fortran_LIBRARIES}
                                    advection_interface_${buildtypedp}  
                                    advection_interface_${buildtypesp} ) 
set_target_properties(${execname_single} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET ${exename_single} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/lib)
