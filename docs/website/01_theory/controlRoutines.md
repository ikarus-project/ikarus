# Control routines

## Load control
A load control object is constructed as follows:
```cpp
auto lc = Ikarus::LoadControl(nonlinearSolver, numLoadSteps, {loadFactorStartValue, loadFactorEndValue});
```
``nonlinearSolver`` is a nonlinear Solver, e.g. Newton-Raphson. ``numLoadSteps`` is the number of load steps, 
``loadFactorStartValue`` is the value of the load factor at the beginning of the simulation (usually 0) and 
``loadFactorEndValue`` is the load factor at the end of the simulation.

The load control is started with the ``run()`` method, i.e. in the example above:
```cpp
lc.run();
```

## Obtaining infos from control routines
The load control is an observable object, i.e. you can subscribe to the messages of the load control.
[Read this page to learn more about the implementation of observer pattern in Ikarus.](observer.md)
The following messages are available:
```cpp
enum class ControlMessages { 
  BEGIN,
  CONTROL_STARTED,
  CONTROL_ENDED,
  STEP_STARTED,
  STEP_ENDED,
  SOLUTION_CHANGED,
  END };
```
