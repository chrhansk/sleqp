# Compiler options
include(CheckCCompilerFlag)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_c_compiler_flag("-fmacro-prefix-map=.=." HAVE_MACRO_PREFIX_MAP)
check_c_compiler_flag("-ffile-prefix-map=.=." HAVE_FILE_PREFIX_MAP)

check_c_source_compiles(
    "
        __attribute__((warn_unused_result))
        void f() {};
        int main(void) {return 0;}
    "
    SLEQP_HAVE_ATTRIBUTE_WARN_UNUSED_RESULT
)

check_c_source_compiles(
    "
        struct [[nodiscard(\"message\")]] error_info { int status; };
        int main(void) {return 0;}
    "
    SLEQP_HAVE_ATTRIBUTE_NODISCARD
    FAIL_REGEX "-Wattributes"
)


check_c_source_compiles(
    "
        __attribute__((__format__(__printf__, 1, 2)))
        void f(const char* format, ...) {};
        int main(void) {return 0;}
    "
    SLEQP_HAVE_ATTRIBUTE_FORMAT
)
