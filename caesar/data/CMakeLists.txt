file(GLOB DATA_FILES 
	${CMAKE_CURRENT_SOURCE_DIR}/*.fits
	${CMAKE_CURRENT_SOURCE_DIR}/*.root
)

set (DATA_FILES ${DATA_FILES} PARENT_SCOPE)
