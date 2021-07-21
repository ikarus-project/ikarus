# Theoretical Background and Implementation Details

This section explains the theoretical background and implementation details 
of various parts of the code. It is dedicated to users who want to extend
or modifiy the implemented functionality or who want to learn more about 
the underlying thoughts and theoretical aspects.

If you are rather interested in the practical use of several aspects, visit
the [tutorial section](../tutorials/tutorialsOverview.md).

## Seperation of interface and implementation
On many of the theory pages you will find a description of an interface and 
a discussion of the implementation. What interface and implementation means is 
explained here with the example of a car.

### Interface of a car
Let's first define the interface of a car. A car is from a certain brand and it 
has a maximum velocity. The interface of a car is then (written as a C++20 
concept):

```cpp
  template <typename CarType>
  concept Car = requires(CarType c) {
    { c.brand() } -> std::same_as<std::string>;
    { c.maxVelocity() } -> std::same_as<double>;
  };
```
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