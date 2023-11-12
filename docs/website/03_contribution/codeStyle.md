# Code style

This section explains some general implementational ideas used in
various parts of the code. It is dedicated to the users who would like to extend
or modify the implemented functionality and/or would like to learn more about
the implementation strategies and certain theoretical aspects.

## General remarks

* The directories and filenames use `camelCase`.
* The source files and the header files have a `cpp` and a `hh` extension, respectively.
* A `clang-format` file is used, which needs to be executed in each extended or modified file before a PR can be merged.
* Readability and value semantics are the essence of the code.
* Classnames use `PascalCase`
* Commenting within the code and the other code styles were influenced by the books by Robert C. Martin[@martinclean] and John K. Ousterhout[@ousterhoutPhilosophySoftwareDesign2021].

\bibliography
