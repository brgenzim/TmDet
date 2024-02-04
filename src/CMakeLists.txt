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

LIST (APPEND FULL_SRC
     ${LibrarySources}
     ${UtilsSources}
     ${TypesSources}
)

ADD_LIBRARY( TmDetLib STATIC ${FULL_SRC} )
INSTALL( TARGETS TmDetLib ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX_LIB} )

ADD_EXECUTABLE( test cli/argTest.cpp)
TARGET_LINK_LIBRARIES(test TmDetLib z)
INSTALL( TARGETS test RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( test2 cli/test2.cpp)
TARGET_LINK_LIBRARIES(test2 z)
INSTALL( TARGETS test2 RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( fragment_cif cli/fragment_cif.cpp)
TARGET_LINK_LIBRARIES(fragment_cif TmDetLib z)
INSTALL( TARGETS fragment_cif RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )


ADD_EXECUTABLE( dssp cli/dssp.cpp)
TARGET_LINK_LIBRARIES(dssp PRIVATE TmDetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS dssp RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( domdist cli/domdist.cpp)
TARGET_LINK_LIBRARIES(domdist PRIVATE z gemmi::gemmi_cpp)
INSTALL( TARGETS domdist RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( domdist_af cli/domdist_af.cpp)
TARGET_LINK_LIBRARIES(domdist_af PRIVATE z gemmi::gemmi_cpp)
INSTALL( TARGETS domdist_af RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

#INSTALL( FILES PdbLib.h DESTINATION include)
