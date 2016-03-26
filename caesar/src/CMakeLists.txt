add_subdirectory(img)
add_subdirectory(imgio)
add_subdirectory(imgfilt)
add_subdirectory(imgproc)
add_subdirectory(utils)

set (SRCS ${SRCS} PARENT_SCOPE)
set (HEADERS ${HEADERS} PARENT_SCOPE)
set (HEADERS_WPATH ${HEADERS_WPATH} PARENT_SCOPE)
set (DICT_SRCS ${DICT_SRCS} PARENT_SCOPE)
set (DICT_HEADERS ${DICT_HEADERS} PARENT_SCOPE)
set (INCLUDE_DIRS ${INCLUDE_DIRS} PARENT_SCOPE)

