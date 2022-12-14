<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

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
* Commenting within the code and the other code styles were influenced by the books by Robert C. Martin[@martinclean] and John K. Ousterhout[@ousterhoutPhilosophySoftwareDesign2021].

## Separation of interface and implementation
The notion of a description of an interface and a discussion of its implementation 
are found in many pages. The question of what an interface and an implementation mean is
answered here with an example of a car.

### Interface of a car
Let's first define the interface of a car. A car is from a certain brand, and it
has a maximum velocity. The car's interface is thus:

- `brand()`: a function which returns the brand name as a `string`
- `maxvelocity()`: a function which returns the maximum velocity as a `double`

This can be written in a more formalized way, e.g., as a C++20 concept, but we currently write it
in this documentation as shown above.

Every object that is then to be declared as a car has to have a member function `brand()` which returns a string, 
and a member function `maxVeloctiy()` which returns a double.

### Implementation of a car
Let's now implement a car.
```cpp
class MyCar
{
  public:
  std::string brand() {return "MyBrand";}
  double maxVelocity() 
  {
    double velocity;
    // calculate maximum velocity with some complicated calculations
    return velocity;
  }
};
```
The class `MyCar` fulfills the car interface and is therefore considered a car.
There can be many other classes that fulfill the car interface. For example, someone could also
implement a class `AnotherCar`.

Thus, the following definitions are summarized:
- Interface: Defines a set of requirements
- Implementation: A specific class which fulfills the interface

\bibliography 