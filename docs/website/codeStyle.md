<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
SPDX-License-Identifier: CC-BY-SA-4.0
-->

### Code style
This section explains the theoretical background and implementation details
of various parts of the code. It is dedicated to users who want to extend
or modify the implemented functionality or who want to learn more about
the implementation thoughts and theoretical aspects.

* The directories and filenames use `camelCase`
* The source files have a `cpp` extension and the header files `hh`.
* We also use a `clang-format` file which needs to be executed in each file before a PR can be merged
* For comments in the code we follow:
!!! cite "Quote of Robert Martin in [VideoLink](https://www.youtube.com/watch?v=2a_ytyt9sf8&list=PLmmYSbUCWJ4x1GO839azG_BBw8rkh-zOj)"
          "The proper use of comments is: 
          To compensate for our failure to express ourselves in code."
       Thus, lots of comments should not be interpreted as good programming style but they should be sign as a failure to express ourselves. 


### Programming style


## Separation of interface and implementation
On many of the theory pages you will find a description of an interface and
a discussion of the implementation. What interface and implementation means is
explained here with the example of a car.

### Interface of a car
Let's first define the interface of a car. A car is from a certain brand and it
has a maximum velocity. The interface of a car is then:

- `brand()`: a function which returns the brand as a string
- `maxvelocity()`: a function which returns the maximum velocity as a double

This can be written in a more formalized way, e.g. as a C++20 concept, but we currently write it
in this documentation as shown above.

Everything that wants to be a car has to have a member function `brand()` which
returns a string and a member function `maxVeloctiy()` which returns a double.

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
The class `MyCar` fulfills the car interface and is a therefore considered a car.
There can be many classes that fulfill the car interface, e.g. someone else could
arrive and implement `class AnotherCar`.

### Summary
- Interface: Defines a set of requirements
- Implementation: A specific class which fulfills the interface

## Member functions and free functions
In this documentation, we also list free functions as a part of the interface. This is indicated
by the arguments in the list of interface functions. An example: If there is something like

- `brand(car)`

in the interface list, this means that there has to be a function which gets a car object as argument
and returns the name of the brand. An implementation for `MyCar` then looks like this:
```cpp
std::string brand(MyCar carObject) {return "MyBrand";}
```
This function is called a free function because it isn't part of the class MyCar but it is free
(and could be defined in another file then the class MyCar).



\bibliography 