
add_subdirectory(sourcefinder)
add_subdirectory(bkgfinder)
add_subdirectory(fits2root)
add_subdirectory(root2fits)
add_subdirectory(saliencyfinder)
add_subdirectory(sourceResidualFinder)
add_subdirectory(sourcefinderMPI)

set (INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)
set (Task_FindSource_SRCS ${Task_FindSource_SRCS} PARENT_SCOPE)
set (Task_FindSourceResidual_SRCS ${Task_FindSourceResidual_SRCS} PARENT_SCOPE)
set (Task_FindSaliency_SRCS ${Task_FindSaliency_SRCS} PARENT_SCOPE)
set (Task_FindBkg_SRCS ${Task_FindBkg_SRCS} PARENT_SCOPE)
set (Task_FITS2ROOT_SRCS ${Task_FITS2ROOT_SRCS} PARENT_SCOPE)
set (Task_ROOT2FITS_SRCS ${Task_ROOT2FITS_SRCS} PARENT_SCOPE)

set (Task_FindSourceMPI_SRCS ${Task_FindSourceMPI_SRCS} PARENT_SCOPE)
