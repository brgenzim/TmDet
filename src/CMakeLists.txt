#/*
# *  CMakeLists.txt
# *  Thu Nov 04 11:12:42 2023
# *  Copyright (c)  2001-2024 Gabor E. Tusnady
# *  <tusi@enzim.hu>
# *  All rights reserved.
# */

AUX_SOURCE_DIRECTORY( . LibrarySources )
AUX_SOURCE_DIRECTORY( cli CliSources )
AUX_SOURCE_DIRECTORY( DTOs DtosSources )
AUX_SOURCE_DIRECTORY( Exceptions ExceptionsSources )
AUX_SOURCE_DIRECTORY( Services ServicesSources )
AUX_SOURCE_DIRECTORY( System SystemSources )
AUX_SOURCE_DIRECTORY( Types TypesSources )
AUX_SOURCE_DIRECTORY( Utils UtilsSources )
AUX_SOURCE_DIRECTORY( ValueObjects ValueObjectsSources )
AUX_SOURCE_DIRECTORY( "${PROJECT_SOURCE_DIR}/../contrib/pugixml/src/" PugiSources)

LIST (APPEND FULL_SRC
     ${LibrarySources}
     ${DtosSources}
     ${ExceptionsSources}
     ${ServicesSources}
     ${SystemSources}
     ${TypesSources}
     ${UtilsSources}
     ${ValueObjectsSources}
     ${PugiSources}
)

ADD_LIBRARY( TmdetLib STATIC ${FULL_SRC} )
INSTALL( TARGETS TmdetLib ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX_LIB} )

ADD_EXECUTABLE( tmdet cli/tmdet.cpp)
TARGET_LINK_LIBRARIES(tmdet PRIVATE TmdetLib z Eigen3::Eigen gemmi::gemmi_cpp curl)
INSTALL( TARGETS tmdet RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )
