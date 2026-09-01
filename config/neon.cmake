message(STATUS "Testing whether ARM NEON code can be used")
set (HAVE_NEON OFF)

try_run(neon_runs neon_compiles
            ${PROJECT_BINARY_DIR}/config
            ${PROJECT_SOURCE_DIR}/config/neon.c)
if(neon_compiles)
    if (neon_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether ARM NEON code can be used -- No (compiles but does not run)")
    else()
        message(STATUS "Testing whether ARM NEON code can be used -- Yes")
        set (HAVE_NEON ON)
    endif()
else()
    try_run(neon_runs neon_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/neon.c
        COMPILE_DEFINITIONS -mfpu=neon)
    if(neon_compiles)
        if (neon_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether ARM NEON code can be used -- No (compiles with -mfpu=neon but does not run)")
        else()
            message(STATUS "Testing whether ARM NEON code can be used -- Yes, with -mfpu=neon")
            add_compile_options(-mfpu=neon)
            set (HAVE_NEON ON)
        endif()
    else()
        message(STATUS "Testing whether ARM NEON code can be used -- No (cannot compile)")
    endif()
endif()