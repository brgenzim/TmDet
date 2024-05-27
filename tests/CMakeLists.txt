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

