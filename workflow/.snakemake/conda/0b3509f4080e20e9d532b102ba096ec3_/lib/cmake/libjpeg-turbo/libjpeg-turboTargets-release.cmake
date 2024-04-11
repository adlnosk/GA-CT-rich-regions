#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libjpeg-turbo::jpeg" for configuration "Release"
set_property(TARGET libjpeg-turbo::jpeg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::jpeg PROPERTIES
  IMPORTED_LOCATION_RELEASE "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libjpeg.so.8.3.2"
  IMPORTED_SONAME_RELEASE "libjpeg.so.8"
  )

list(APPEND _cmake_import_check_targets libjpeg-turbo::jpeg )
list(APPEND _cmake_import_check_files_for_libjpeg-turbo::jpeg "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libjpeg.so.8.3.2" )

# Import target "libjpeg-turbo::turbojpeg" for configuration "Release"
set_property(TARGET libjpeg-turbo::turbojpeg APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::turbojpeg PROPERTIES
  IMPORTED_LOCATION_RELEASE "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libturbojpeg.so.0.3.0"
  IMPORTED_SONAME_RELEASE "libturbojpeg.so.0"
  )

list(APPEND _cmake_import_check_targets libjpeg-turbo::turbojpeg )
list(APPEND _cmake_import_check_files_for_libjpeg-turbo::turbojpeg "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libturbojpeg.so.0.3.0" )

# Import target "libjpeg-turbo::turbojpeg-static" for configuration "Release"
set_property(TARGET libjpeg-turbo::turbojpeg-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::turbojpeg-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "ASM_NASM;C"
  IMPORTED_LOCATION_RELEASE "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libturbojpeg.a"
  )

list(APPEND _cmake_import_check_targets libjpeg-turbo::turbojpeg-static )
list(APPEND _cmake_import_check_files_for_libjpeg-turbo::turbojpeg-static "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libturbojpeg.a" )

# Import target "libjpeg-turbo::jpeg-static" for configuration "Release"
set_property(TARGET libjpeg-turbo::jpeg-static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libjpeg-turbo::jpeg-static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "ASM_NASM;C"
  IMPORTED_LOCATION_RELEASE "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libjpeg.a"
  )

list(APPEND _cmake_import_check_targets libjpeg-turbo::jpeg-static )
list(APPEND _cmake_import_check_files_for_libjpeg-turbo::jpeg-static "/work/project/briefwp3/Adela/GA_CT/workflow/.snakemake/conda/0b3509f4080e20e9d532b102ba096ec3_/lib/libjpeg.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
