cmake_minimum_required (VERSION 3.14.0)
project(advectionHPCtester Fortran )
include( CMakeParseArguments )
enable_language (Fortran)

set(CMAKE_Fortran_FLAGS "-g ")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(optsprecsp -DFLOAT_PRECISION=4) 
set(optsprecdp -DFLOAT_PRECISION=8) 
set(optswithgpu -DCUDACODE -DTIMERSCUDA )
set(optswithmpinogpu -DNOCUDACODE -DTIMERSMPI )
set(optsnompinogpu -DNOCUDACODE -DTIMERSCPU )

option(ENABLE_GPU "Configure for NVIDIA GPU" ON)

if(ENABLE_GPU)
find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
  message(STATUS "Building with GPU support")
  enable_language (CUDA)
  find_library(LIBNVTOOLSEXT nvToolsExt PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
  find_library(LIBPNETCDF libpnetcdf PATHS "/usr/lib/x86_64-linux-gnu")
  #  set(optsgpu -DCUDACODE -DTIMERSCUDA  )
else()
  message(STATUS "Building without GPU support")
  #  set(optsgpu -DNOCUDACODE   )
endif()

option(ENABLE_MPI "Configure for MPI parallelization" OFF)

if(ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  message(STATUS "Building with MPI")
    set(optsmpi  -DPUREMPI -DMPI90 )
  message (STATUS ${MPI_Fortran_INCLUDE_PATH})
  include_directories(${MPI_Fortran_INCLUDE_PATH})
else()
  message(STATUS "Building without MPI")
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

option(ENABLE_FAST "Optimize - ON, Debug - OFF"  OFF)

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
       set(optims  -O1 -Mpreprocess )
       set(debugs  -C traceback -Mbounds -Mchkfpstk -Kieee -Mlist -Mpreprocess )
  endif()

endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debugs}")
set(CMAKE_Fortran_FLAGS_TMP  "${optims} ${optsmpi} ${optpnet} ${optstest}")
#if(ENABLE_GPU)
#set(CMAKE_Fortran_FLAGS_TMP "${CMAKE_Fortran_FLAGS_TMP} -cuda")
#endif()

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
#mpdata3d-rho-gauge-gpubc.f90
velprd_driver.f90
module_velprd.f90
mpdataop.f90
)

set(advection_drv_sources_preprocess ${advection_drv_sources_preprocess}
advec.f90
advec_init_finalize.f90
#advection_CEtimestep_twostep_GPUBC.f90
)

set(advection_drvint_sources_preprocess ${advection_drvint_sources_preprocess}
advec_interface_sp.f90
advec_init_interface_sp.f90
advec_interface_dp.f90
advec_init_interface_dp.f90
advec_interface_spgpu.f90
advec_init_interface_spgpu.f90
advec_interface_dpgpu.f90
advec_init_interface_dpgpu.f90
)

set(advection_tst_sources_preprocess ${advection_tst_sources_preprocess}
advection_test_setup.f90
)

function(DEFINE_PREPROCESS_SED) #subdir_arg destdir_arg files_arg custopt_arg)
set (options ONGPU )
set (oneValueArgs SUBDIR BUILDTYPE)
set (multiValueArgs SOURCEFILES SPECOPTS)
cmake_parse_arguments(DPSPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(subdir  ${DPSPREFIX_SUBDIR})
set(destdir ${DPSPREFIX_BUILDTYPE})
set(custopt ${DPSPREFIX_SPECOPTS})
set(files   ${DPSPREFIX_SOURCEFILES})
message (STATUS "Subdir: " ${subdir})
message (STATUS "Buildtype: " ${destdir})
message (STATUS "Specopts: " ${custopt})
message (STATUS "Files " ${files})
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/${destdir}/${subdir}") 
set(localdir "${PROJECT_SOURCE_DIR}/${subdir}")
set(destdir  "${destdir}/${subdir}")
foreach( file_i ${files})
        message (STATUS "define preprocessing for file: ${destdir}/${file_i}")
        add_custom_command( 
		COMMAND gcc -E -P ${optsmpi} ${optstest} ${custopt} -cpp  ${localdir}/${file_i} -DSYMBOLSUFFIX=${DPSPREFIX_BUILDTYPE}
        | sed -e "s/\&/\&\\n/g" -e "s/\;/\;\\n/g" -e "s/>>>/>>>\\n/g" -e "/^$/d" > "${destdir}/${file_i}.pp.f90"
        DEPENDS  ${localdir}/${file_i}
        OUTPUT  ${destdir}/${file_i}.pp.f90
        COMMENT "Preprocessing ${localdir}/${file_i} to ${destdir}/${file_i}.pp.f90"
        VERBATIM
)
endforeach( file_i )
endfunction()

function(GENERATE_BUILD) 
set (options ONGPU )
set (oneValueArgs BUILDTYPE)
set (multiValueArgs SPECOPTS DRIVERFILES)
cmake_parse_arguments(GBPREFIX "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
set(buildtype ${GBPREFIX_BUILDTYPE}) #buildtype_arg
set(specopts  ${GBPREFIX_SPECOPTS}) #specopts_arg} )
set(driverfile ${GBPREFIX_DRIVERFILES}) #  ${driverfile_arg} )
set(ongpu ${GBPREFIX_ONGPU}) #  ${driverfile_arg} )
message (STATUS "Buildtype: " ${buildtype})
message (STATUS "Specopts: " ${specopts})
message (STATUS "Driver file: " ${driverfile})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/common) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/common/argsparser) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/advection) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/include) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/lib) 
set(advection_sources )  

define_preprocess_sed(SUBDIR "common"                   BUILDTYPE  ${buildtype}  SOURCEFILES ${common_sources}                   SPECOPTS ${specopts})
define_preprocess_sed(SUBDIR "common"                   BUILDTYPE  ${buildtype}  SOURCEFILES ${common_sources_preprocess}        SPECOPTS ${specopts})
define_preprocess_sed(SUBDIR "advection/src_algorithms" BUILDTYPE  ${buildtype}  SOURCEFILES ${advection_alg_sources_preprocess} SPECOPTS ${specopts})
define_preprocess_sed(SUBDIR "advection/src_drivers"    BUILDTYPE  ${buildtype}  SOURCEFILES ${advection_drv_sources_preprocess} SPECOPTS ${specopts})
define_preprocess_sed(SUBDIR "advection/src_drivers"    BUILDTYPE  ${buildtype}  SOURCEFILES ${driverfile} SPECOPTS ${specopts})
define_preprocess_sed(SUBDIR "advection/src_testdata"   BUILDTYPE  ${buildtype}  SOURCEFILES ${advection_tst_sources_preprocess} SPECOPTS ${specopts})

list(TRANSFORM advection_alg_sources_preprocess PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/advection/src_algorithms/)
list(TRANSFORM advection_alg_sources_preprocess APPEND .pp.f90)
list(TRANSFORM advection_drv_sources_preprocess PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/advection/src_drivers/)
list(TRANSFORM advection_drv_sources_preprocess APPEND .pp.f90)
list(TRANSFORM driverfile PREPEND ${CMAKE_BINARY_DIR}/${buildtype}/advection/src_drivers/)
list(TRANSFORM driverfile APPEND .pp.f90)
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
if(ongpu)
	target_compile_options(   common_${buildtype} PRIVATE -cuda -Wl,--exclude-libs,ALL ) 
	target_compile_options(advection_${buildtype} PRIVATE -cuda -Wl,--exclude-libs,ALL ) 
endif()
target_include_directories(common_${buildtype} PRIVATE ${commondir} ) 
target_include_directories(advection_${buildtype} PRIVATE ${advectiondir}) 

#list(TRANSFORM driverfile PREPEND ${CMAKE_SOURCE_DIR}/advection/src_drivers/)
if(ongpu)
	add_library(advection_interface_${buildtype} STATIC ${driverfile})
else()
	add_library(advection_interface_${buildtype} SHARED ${driverfile})
endif()
add_dependencies(advection_interface_${buildtype}  advection_${buildtype} )
target_include_directories(advection_interface_${buildtype} PRIVATE ${commondir} ) 
target_include_directories(advection_interface_${buildtype} PRIVATE ${advectiondir}) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/include) 
set_target_properties(advection_interface_${buildtype} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set_target_properties(advection_interface_${buildtype} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
if(ongpu)
	target_compile_options(             common_${buildtype} PRIVATE -cuda -gpu=rdc  -Wl,--exclude-libs,ALL) 
	target_compile_options(          advection_${buildtype} PRIVATE -cuda -gpu=rdc  -Wl,--exclude-libs,ALL) 
	target_compile_options(advection_interface_${buildtype} PRIVATE -cuda -gpu=rdc  -Wl,--exclude-libs,ALL) 
        target_link_options(common_${buildtype}  PUBLIC  -lnvToolsExt -r -cuda -Wl,--exclude-libs,ALL )
        target_link_options(advection_${buildtype}  PUBLIC  -lnvToolsExt -r  -cuda -Wl,--exclude-libs,ALL )
        target_link_options(advection_interface_${buildtype}  PUBLIC  -lnvToolsExt -r -cuda -Wl,--exclude-libs,ALL )
	#	target_link_options(advection_${buildtype} PUBLIC -Xcompiler -fvisibility=hidden ) 
	#	target_link_options(common_${buildtype} PUBLIC -Xcompiler -fvisibility=hidden ) 
endif()

endfunction(GENERATE_BUILD)

set(buildtype sp)
set(driverfile 
    advec_interface_${buildtype}.f90;
    advec_init_interface_${buildtype}.f90
    )
set(custopts ${optsprecsp} "${optsnompinogpu}")
generate_build( SPECOPTS  ${custopts} DRIVERFILES ${driverfile} BUILDTYPE ${buildtype} )
set(execname_single advectiontest.${buildtype})
add_link_options( -ldl)
#if(ENABLE_GPU)
#add_link_options( -cuda -lnvToolsExt)
#set(execname_single advectiontest.gpu.${buildtype})
#else()
#set(execname_single advectiontest.cpu.${buildtype})
#endif()
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal_${buildtype}.f90")
add_dependencies(${execname_single} common_${buildtype} advection_${buildtype} )
target_link_libraries(advection_interface_${buildtype} PRIVATE advection_${buildtype})
target_link_libraries(${execname_single}  advection_interface_${buildtype})
target_include_directories(${execname_single} PRIVATE ${CMAKE_BINARY_DIR}/include) 
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)

set(buildtype dp)
set(driverfile 
    advec_interface_${buildtype}.f90;
    advec_init_interface_${buildtype}.f90)
set(custopts ${optsprecdp} "${optsnompinogpu}")
#generate_build(${buildtype} ${custopts}  ${driverfile}) 
#generate_build(BUILDTYPE ${buildtype} SPECOPTS  ${custopts} DRIVERFILES  ${driverfile}) 
generate_build( SPECOPTS  ${custopts} DRIVERFILES ${driverfile} BUILDTYPE ${buildtype} )
set(execname_single advectiontest.${buildtype})
#if(ENABLE_GPU)
#add_link_options( -cuda -lnvToolsExt)
#set(execname_single advectiontest.gpu.${buildtype})
#else()
#set(execname_single advectiontest.cpu.${buildtype})
#endif()

add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal_${buildtype}.f90")
add_dependencies(${execname_single} common_${buildtype} advection_${buildtype} )
target_link_libraries(advection_interface_${buildtype} PRIVATE  advection_${buildtype})
target_link_libraries(${execname_single}  advection_interface_${buildtype})
target_include_directories(${execname_single} PRIVATE ${CMAKE_BINARY_DIR}/include) # ${solverdir} ${diffusiondir})
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)

set(buildtype dpgpu)
set(driverfile 
    advec_interface_${buildtype}.f90;
    advec_init_interface_${buildtype}.f90)
set(custopts ${optsprecdp} "${optswithgpu}")
#generate_build(${buildtype} @custopts@  ${driverfile}) 
#generate_build(BUILDTYPE ${buildtype} SPECOPTS  ${custopts} DRIVERFILES  ${driverfile}) 
generate_build( ONGPU TRUE SPECOPTS  ${custopts} DRIVERFILES ${driverfile} BUILDTYPE ${buildtype} )
set(execname_single advectiontest.${buildtype})

add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal_${buildtype}.f90")
target_compile_options(${execname_single} PRIVATE -cuda -gpu=rdc ${optswithgpu}) 
add_dependencies(${execname_single} common_${buildtype} advection_${buildtype} )
target_link_libraries(advection_interface_${buildtype} PRIVATE  advection_${buildtype})
target_link_options(${execname_single} PRIVATE -cuda )# -lnvToolsExt)
target_link_libraries(${execname_single}  advection_interface_${buildtype})
if(ENABLE_GPU)
	target_compile_options(${execname_single} PRIVATE -cuda) 
endif()
target_include_directories(${execname_single} PRIVATE ${CMAKE_BINARY_DIR}/include) # ${solverdir} ${diffusiondir})
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)

set(buildtype spgpu)
set(driverfile 
    advec_interface_${buildtype}.f90;
    advec_init_interface_${buildtype}.f90)
set(custopts ${optsprecsp} "${optswithgpu}")
generate_build( ONGPU TRUE SPECOPTS  ${custopts} DRIVERFILES ${driverfile} BUILDTYPE ${buildtype} )
set(execname_single advectiontest.${buildtype})
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal_${buildtype}.f90")
target_compile_options(${execname_single} PRIVATE -cuda -gpu=rdc ${optswithgpu}) 
add_dependencies(${execname_single} common_${buildtype} advection_${buildtype} )
target_link_options(${execname_single} PRIVATE -cuda)# -lnvToolsExt)
target_link_libraries(advection_interface_${buildtype} PRIVATE  advection_${buildtype})
target_link_libraries(${execname_single}  advection_interface_${buildtype})
if(ENABLE_GPU)
	target_compile_options(${execname_single} PRIVATE -cuda) 
endif()
target_include_directories(${execname_single} PRIVATE ${CMAKE_BINARY_DIR}/include) # ${solverdir} ${diffusiondir})
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)


set(buildtype uni)
set(buildtypesp sp)
set(buildtypedp dp)
set(buildtypespgpu spgpu)
set(buildtypedpgpu dpgpu)
set(execname_single advectiontest.uni)
add_executable(${execname_single} "advection/src_testapps/advection_dwarf_test_minimal.f90")
target_compile_options(${execname_single} PRIVATE -r -cuda -gpu=rdc ${optswithgpu}) 
add_dependencies(${execname_single} advection_interface_${buildtypedp} 
                                    advection_interface_${buildtypesp}   
                                    advection_interface_${buildtypedpgpu}   
                                    advection_interface_${buildtypespgpu}  )
			    #target_link_options(${execname_single} PUBLIC  -lnvToolsExt  -cuda -Wl,--exclude-libs,ALL )
target_link_libraries(${execname_single} PRIVATE ${MPI_Fortran_LIBRARIES}
	#                                    advection_interface_${buildtypedp}  
	#                            advection_interface_${buildtypesp}
	                           advection_interface_${buildtypedpgpu}  
                                   advection_interface_${buildtypespgpu}
                                    )
if(ENABLE_GPU)
	#	target_compile_options(${execname_single} PRIVATE -cuda -gpu=rdc -Wl,--exclude-libs,ALL) 
endif()
set_target_properties(${execname_single} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set_property(TARGET ${exename_single} PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET ${exename_single} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${buildtype}/lib)
