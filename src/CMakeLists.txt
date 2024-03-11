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
AUX_SOURCE_DIRECTORY( Types TypesSources )
AUX_SOURCE_DIRECTORY( TmdetLib TmdetLibSources )
AUX_SOURCE_DIRECTORY( "${PROJECT_SOURCE_DIR}/../contrib/pugixml/src/" PugiSources)

LIST (APPEND FULL_SRC
     ${LibrarySources}
     ${UtilsSources}
     ${TypesSources}
     ${TmdetLibSources}
     ${PugiSources}
)

ADD_LIBRARY( TmdetLib STATIC ${FULL_SRC} )
INSTALL( TARGETS TmdetLib ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX_LIB} )

#ADD_EXECUTABLE( test cli/argTest.cpp)
#TARGET_LINK_LIBRARIES(test TmdetLib z)
#INSTALL( TARGETS test RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( tmdet cli/tmdet.cpp)
TARGET_LINK_LIBRARIES(tmdet PRIVATE TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS tmdet RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

#ADD_EXECUTABLE( fragment_cif cli/fragment_cif.cpp)
#TARGET_LINK_LIBRARIES(fragment_cif TmdetLib z)
#INSTALL( TARGETS fragment_cif RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( dssp cli/dssp.cpp)
TARGET_LINK_LIBRARIES(dssp PRIVATE TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS dssp RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

#ADD_EXECUTABLE( domdist cli/domdist.cpp)
#TARGET_LINK_LIBRARIES(domdist PRIVATE z gemmi::gemmi_cpp)
#INSTALL( TARGETS domdist RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

#ADD_EXECUTABLE( domdist_af cli/domdist_af.cpp)
#TARGET_LINK_LIBRARIES(domdist_af PRIVATE z gemmi::gemmi_cpp)
#INSTALL( TARGETS domdist_af RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

#INSTALL( FILES PdbLib.h DESTINATION include)
