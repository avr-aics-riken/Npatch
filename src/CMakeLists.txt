
include_directories("${PROJECT_BINARY_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(Npatch Npt.cxx FNpt.cxx)

add_definitions("${REAL_OPT}")

########### install files ###############

install(TARGETS Npatch DESTINATION  ${PROJECT_NAME}/lib)

install(FILES ../include/CalcGeo.h ../include/CalcGeo_Matrix.h
              ../include/FNpt.h ../include/Npt.h 
              ${PROJECT_BINARY_DIR}/include/npt_Version.h
	      DESTINATION ${PROJECT_NAME}/include)
