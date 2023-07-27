# Extract value of a '"#define X 1"' from a header file
function(define_extract_int file name)
  file(STRINGS "${file}"
    file_result
    REGEX "^#define ${name} +[0-9]+")
  string(REGEX REPLACE "^#define ${name} +([0-9]+)" "\\1" replace_result ${file_result})
  set(${name} ${replace_result} PARENT_SCOPE)
endfunction()

# Extract value of a '#define X "1.1.1"' from a header file
function(define_extract_version file name)
  file(STRINGS "${file}"
    file_result
    REGEX "^#define ${name} \"[0-9\.]+\"")
  string(REGEX REPLACE "^#define ${name} \"([0-9\.]+)\"" "\\1" replace_result ${file_result})
  set(${name} ${replace_result} PARENT_SCOPE)
endfunction()
