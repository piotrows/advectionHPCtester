cmake_minimum_required (VERSION 3.14.0)
project(advectionHPCtester Fortran )
enable_language (Fortran)
set(CMAKE_Fortran_FLAGS "-g")
#execute_process(COMMAND bash  "${PROJECT_SOURCE_DIR}/advection/src_algorithms/prepgpu.sh" 
#  	        WORKING_DIRECTORY  "${PROJECT_SOURCE_DIR}/advection/src_algorithms/" )
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

option(ENABLE_MPI "Configure for MPI parallelization" ON)

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
    set(optims  -Ofast -cpp )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(optims  -O3 -xHost -unroll-aggressive -qopt-prefetch -align commons )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(optims  -O2 -Mpreprocess )
  endif()

else()

  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
       set(optims  -ffree-form -fimplicit-none -cpp  )
       set(debugs  -g -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow  -finit-real=snan  )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
       set(optims  -O1 )
       set(debugs  -g -traceback -check bounds -ftrapuv -init=snan -implicitnone -fpe0 -ftz )
  endif()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
       set(optims  -O1 -Mpreprocess )
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

set(commondir    "${PROJECT_SOURCE_DIR}/include/common")
set(advectiondir "${PROJECT_SOURCE_DIR}/include/advection")

set(common_sources_to_preprocess ${common_sources_to_preprocess}
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

set(advection_alg_sources_to_preprocess ${advection_alg_sources_to_preprocess}
upwind3d-gpubc.f90
antidiff3d-gauge-gpubc.f90
antidiff3d-standard-gpubc.f90
mpdata3d-rho-gauge-gpubc.f90
velprd_driver.f90
module_velprd.f90
)

set(advection_drv_sources_to_preprocess ${advection_drv_sources_to_preprocess}
advec.f90
advection_CEtimestep_twostep_GPUBC.f90
)
set(advection_tst_sources_to_preprocess ${advection_tst_sources_to_preprocess}
advection_test_setup.f90
)


##foreach( file_i ${common_sources_to_preprocess})
## execute_process(COMMAND cp ${file_i} ${file_i}.pp.f90
##               WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/common )
##endforeach( file_i )
#foreach( file_i ${advection_alg_sources_to_preprocess})
# execute_process(COMMAND cp ${file_i} ${file_i}.pp.f90
#               WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/advection/src_algorithms )
#endforeach( file_i )
#foreach( file_i ${advection_drv_sources_to_preprocess})
# execute_process(COMMAND cp ${file_i} ${file_i}.pp.f90
#               WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/advection/src_drivers )
#endforeach( file_i )
#foreach( file_i ${advection_tst_sources_to_preprocess})
# execute_process(COMMAND cp ${file_i} ${file_i}.pp.f90
#               WORKING_DIRECTORY  ${PROJECT_SOURCE_DIR}/advection/src_testdata )
#endforeach( file_i )

add_library(commontmp INTERFACE  )
add_library(advectiontmp INTERFACE  )
set(localdir "${PROJECT_SOURCE_DIR}/common")
foreach( file_i ${common_sources_to_preprocess})
        message (STATUS "inside loop:${localdir}/${file_i}")
        add_custom_command(
        COMMAND gcc-13 -E -P ${optsgpu} ${optsmpi}  -cpp  ${localdir}/${file_i}
        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${localdir}/${file_i}.pp.f90"
        DEPENDS  ${localdir}/${file_i}
        OUTPUT  ${localdir}/${file_i}.pp.f90
        COMMENT "Preprocessing ${PROJECT_SOURCE_DIR}/common/${file_i} to ${localdir}/${file_i}.pp.f90"
        VERBATIM
)
#        add_custom_target( ${file_i}.target
#        COMMAND gcc-13 -E -P ${optsgpu} ${optsmpi}  -cpp  ${localdir}/${file_i}
#        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${localdir}/${file_i}.pp.f90"
#        DEPENDS  ${localdir}/${file_i}
#        COMMENT "Preprocessing ${PROJECT_SOURCE_DIR}/common/${file_i} to ${localdir}/${file_i}.pp.f90"
#        VERBATIM
#)
#add_dependencies (commontmp ${file_i}.target)
endforeach( file_i )

set(localdir "${PROJECT_SOURCE_DIR}/advection/src_algorithms")
foreach( file_i ${advection_alg_sources_to_preprocess})
        message (STATUS "inside loop:${localdir}/${file_i}")
        add_custom_target( ${file_i}.target
        COMMAND gcc-13 -E -P ${optsgpu} ${optsmpi} -cpp  ${localdir}/${file_i}
        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${localdir}/${file_i}.pp.f90"
        DEPENDS  ${localdir}/${file_i}
        COMMENT "APreprocessing ${PROJECT_SOURCE_DIR}/common/${file_i} to ${localdir}/${file_i}.pp.f90"
        VERBATIM
)
add_dependencies (advectiontmp ${file_i}.target)
endforeach( file_i )

set(localdir "${PROJECT_SOURCE_DIR}/advection/src_drivers")
foreach( file_i ${advection_drv_sources_to_preprocess})
        message (STATUS "inside loop:${localdir}/${file_i}")
        add_custom_target( ${file_i}.target
        COMMAND gcc-13 -E -P ${optsgpu} ${optsmpi} -cpp  ${localdir}/${file_i}
        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${localdir}/${file_i}.pp.f90"
        DEPENDS  ${localdir}/${file_i}
        COMMENT "BPreprocessing ${PROJECT_SOURCE_DIR}/common/${file_i} to ${localdir}/${file_i}.pp.f90"
        VERBATIM
)
add_dependencies (advectiontmp ${file_i}.target)
endforeach( file_i )

set(localdir "${PROJECT_SOURCE_DIR}/advection/src_testdata")
foreach( file_i ${advection_tst_sources_to_preprocess})
        message (STATUS "inside loop:${localdir}/${file_i}")
        add_custom_target( ${file_i}.target
        COMMAND gcc-13 -E -P ${optsmpi} ${optsgpu} -cpp  ${localdir}/${file_i}
        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${localdir}/${file_i}.pp.f90"
        DEPENDS  ${localdir}/${file_i}
        COMMENT "CPreprocessing ${PROJECT_SOURCE_DIR}/common/${file_i} to ${localdir}/${file_i}.pp.f90"
        VERBATIM
)
add_dependencies (advectiontmp ${file_i}.target)
endforeach( file_i )


set(advection_sources ${advection_sources})
list(TRANSFORM advection_alg_sources_to_preprocess PREPEND advection/src_algorithms/)
list(TRANSFORM advection_alg_sources_to_preprocess APPEND .pp.f90)
list(TRANSFORM advection_drv_sources_to_preprocess PREPEND advection/src_drivers/)
list(TRANSFORM advection_drv_sources_to_preprocess APPEND .pp.f90)
list(TRANSFORM advection_tst_sources_to_preprocess PREPEND advection/src_testdata/)
list(TRANSFORM advection_tst_sources_to_preprocess APPEND .pp.f90)
list(APPEND advection_sources ${advection_alg_sources_to_preprocess} )
list(APPEND advection_sources ${advection_drv_sources_to_preprocess} )
list(APPEND advection_sources ${advection_tst_sources_to_preprocess} )
set(common_sources ${common_sources}
precisions.f90
argsparser
bconditions.f90
epsilons.f90
eulag_lib_parameters.f90
noise.f90
argsparser/argsparser.f90)

list(TRANSFORM common_sources PREPEND common/)
list(TRANSFORM common_sources_to_preprocess PREPEND common/)
list(TRANSFORM common_sources_to_preprocess APPEND .pp.f90)
list(APPEND common_sources ${common_sources_to_preprocess} )

add_library(common    STATIC ${common_sources})   
add_library(advection STATIC ${advection_sources})
add_dependencies(common commontmp )
add_dependencies(advection common advectiontmp commontmp )

set_target_properties(common    PROPERTIES Fortran_MODULE_DIRECTORY ${commondir})
set_target_properties(advection PROPERTIES Fortran_MODULE_DIRECTORY ${advectiondir})
include_directories(${commondir} ${advectiondir}) # ${solverdir} ${diffusiondir})

link_libraries(${MPI_Fortran_LIBRARIES} common advection)
if(ENABLE_GPU)
add_link_options( -cuda -lnvToolsExt)
set(execname_single advectiontest.gpu)
#set(execname_fullset tladvectiontest.gpu)
else()
set(execname_single advectiontest.cpu)
set(execname_fullset tladvectiontest.cpu)
endif()
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_modular.f90")
add_dependencies(${execname_single} common advection )
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)
