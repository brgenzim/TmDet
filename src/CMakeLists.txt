#/*
# *  CMakeLists.txt
# *  Thu Nov 04 11:12:42 2023
# *  Copyright (c)  2001-2024 Gabor E. Tusnady
# *  <tusi@enzim.hu>
# *  All rights reserved.
# */

AUX_SOURCE_DIRECTORY( . LibrarySources )
AUX_SOURCE_DIRECTORY( cli CliSources )
AUX_SOURCE_DIRECTORY( Utils UtilsSources )
AUX_SOURCE_DIRECTORY( DTOs DtoSources )
AUX_SOURCE_DIRECTORY( Services ServicesSources )
AUX_SOURCE_DIRECTORY( Types TypeSources )
AUX_SOURCE_DIRECTORY( "${PROJECT_SOURCE_DIR}/../contrib/pugixml/src/" PugiSources)

LIST (APPEND FULL_SRC
     ${LibrarySources}
     ${UtilsSources}
     ${DtoSources}
     ${PugiSources}
     ${ServicesSources}
     ${TypeSources}
)

ADD_LIBRARY( TmdetLib STATIC ${FULL_SRC} )
INSTALL( TARGETS TmdetLib ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX_LIB} )

#ADD_EXECUTABLE( test cli/argTest.cpp)
#TARGET_LINK_LIBRARIES(test TmdetLib z)
#INSTALL( TARGETS test RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( tmdet cli/tmdet.cpp)
TARGET_LINK_LIBRARIES(tmdet PRIVATE TmdetLib z gemmi::gemmi_cpp Eigen3::Eigen)
INSTALL( TARGETS tmdet RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )


ADD_EXECUTABLE( fragment_cif cli/fragment_cif.cpp)
TARGET_LINK_LIBRARIES(fragment_cif TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS fragment_cif RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

#INSTALL( FILES PdbLib.h DESTINATION include)
