# useful macros/fcts
# imported from waves project
# Authors Romain Boman


MACRO(MACRO_AddTest srcDir)
    #message(STATUS "Adding test directory ${srcDir}")
    file(GLOB tfiles RELATIVE ${srcDir} ${srcDir}/*)
    #message(STATUS "tfiles=${tfiles}")
    foreach(tfile ${tfiles})
        set(spath ${srcDir}/${tfile})
        if((NOT IS_DIRECTORY ${spath}) AND (${spath} MATCHES ".+\\fsi.py$" OR ${spath} MATCHES ".+\\fsi_adj.py$"))
            string(REPLACE "${PROJECT_SOURCE_DIR}/" "" strip ${spath}) 
            message(STATUS "Adding test ${strip}")
            add_test(NAME ${strip} 
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} 
                     COMMAND ${PYTHON_EXECUTABLE} run.py ${strip})
        else()
            MACRO_AddTest(${srcDir}/${tfile})
        endif()
    endforeach()
ENDMACRO()
