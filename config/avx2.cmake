message(STATUS "Testing whether avx2 code can be used")
set (HAVE_AVX2 OFF)

try_run(avx2_runs avx2_compiles
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/avx2.c)
if(avx2_compiles)
    if (avx2_runs MATCHES FAILED_TO_RUN)
        message(STATUS "Testing whether avx2 code can be used -- No (compiles but does not run)")
    else()
        message(STATUS "Testing whether avx2 code can be used -- Yes")
        set (HAVE_AVX2 ON)
    endif()
else()

    # 2nd try
    try_run(avx2_runs avx2_compiles
        ${PROJECT_BINARY_DIR}/config
        ${PROJECT_SOURCE_DIR}/config/avx2.c
        COMPILE_DEFINITIONS -mavx2)
    if(avx2_compiles)
        if(avx2_runs MATCHES FAILED_TO_RUN)
            message(STATUS "Testing whether avx2 code can be used -- No (compiles with -mavx2 but does not run)")
        else()
            message(STATUS "Testing whether avx2 code can be used -- Yes, with -mavx2")
            add_compile_options(-mavx2)
            set (HAVE_AVX2 ON)
        endif()
    else()
        message(STATUS "Testing whether avx2 code can be used -- No")
    endif()
endif()