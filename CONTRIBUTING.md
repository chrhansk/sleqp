# Contributing

To contribute documentation or code, please send a pull request.
Please ensure that any code

* compiles correctly (the CI pipeline passes)
* is sufficiently tested (does not negatively impact coverage)
* is formatted correctly (see below)

## Code formatting

We use [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
to enforce a consistent style, according to the `.clang-format` file
in the root directory. To format your code, run

```
clang-format -style=file path/to/source.c
```

in any project folder.

## Licensing

This project is distributed according to the `LICENSE`
file. Consequently, we expect all contributors to agree that
contributions will eventually become part of the project subjecting to
this license as well.
