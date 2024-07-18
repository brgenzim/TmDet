#/*
# *  CMakeLists.txt
# *  Mon May 27 22:17:57 2024
# *  Copyright (c)  2001-2024 Gabor E. Tusnady
# *  <tusi@enzim.hu>
# *  All rights reserved.
# */

AUX_SOURCE_DIRECTORY( . LibrarySources )
AUX_SOURCE_DIRECTORY( units UnitsSources )
AUX_SOURCE_DIRECTORY( features FeaturesSources )

LIST (APPEND FULL_SRC
     ${UnitsSources}
     ${FeaturesSources}
)

ADD_EXECUTABLE( oligomerTest units/oligomerTest.cpp)
TARGET_LINK_LIBRARIES( oligomerTest PRIVATE TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS oligomerTest RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( secondaryStructVectorTest units/secondaryStructVectorTest.cpp)
TARGET_LINK_LIBRARIES( secondaryStructVectorTest PRIVATE TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS secondaryStructVectorTest RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( alignmentTest units/alignmentTest.cpp)
TARGET_LINK_LIBRARIES( alignmentTest PRIVATE TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS alignmentTest RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )

ADD_EXECUTABLE( dsspTest units/dsspTest.cpp)
TARGET_LINK_LIBRARIES( dsspTest PRIVATE TmdetLib z gemmi::gemmi_cpp)
INSTALL( TARGETS dsspTest RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX_BIN} )
