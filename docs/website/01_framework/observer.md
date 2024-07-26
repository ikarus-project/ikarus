---
status: new
---

# Observer and Observable

To write output messages when desired by the user, an observer pattern is implemented in Ikarus.
Four things are necessary to understand the implementation of observer patterns: `Messages`, `IObservable`,
`IObserver` and `Subscriptions`.

## Messages

A message class is a list of possible events that can happen and might be of interest. The messages that are used for
nonlinear solvers are listed below as an example.

```cpp
enum class NonLinearSolverMessages {
  BEGIN,
  INIT,
  ITERATION_STARTED,
  ITERATION_ENDED,
  RESIDUALNORM_UPDATED,
  CORRECTIONNORM_UPDATED,
  SOLUTION_CHANGED,
  SOLVER_FINISHED,
  END
};
```

## IObservable

A class can be observable.
The class then sends notifications when certain events are happening.
To become observable, a class
must inherit from `IObservable<MessageType, StateType>`.
Here, `MessageType`
is the `enum` of messages to use (see above)
and `StateType` can either be `NonLinearSolverState` or `ControlState`.

There are two aliases available for the class `IObservable`:

```cpp
using ControlObservable = IObservable<ControlMessages, ControlState>;
using NonLinearSolverObservable = IObservable<NonLinearSolverMessages, NonLinearSolverState>;
```

`IObservable` must adhere to the following concept:

```cpp
template <typename MT, typename ST>
concept Observable = ObserverMessage<MT> and ObserverState<ST>;
```

with

```cpp
template <typename MT>
concept ObserverMessage = std::is_same_v<MT, ControlMessages> or std::is_same_v<MT, NonLinearSolverMessages>;

template <typename ST>
concept ObserverState = std::is_same_v<ST, ControlState> or std::is_same_v<ST, NonLinearSolverState>;
```

## IObserver

A class can be an observer.
The class is then notified when events are happening and can perform actions.
To become an observer, the class must inherit from `IObserver<IObservable<MessageType, StateType>>`.
A basic example is shown below.

```cpp
class OurFirstObserver : public IObserver<NonLinearSolverObservable> {
public:
  using MessageType = typename NonLinearSolverObservable::MessageType;
  using StateType   = typename NonLinearSolverObservable::StateType;
  void updateImpl(MessageType message, const StateType&) override {
    if (message == NonLinearSolverMessages::ITERATION_STARTED) std::cout << "Iteration started.\n";
  }
};
```

The observer has to implement the function `void updateImpl(MessageType message, const StateType& state)`.
In this function, all actions can
be implemented that should be performed when the corresponding message is received.

## Subscriptions

There are a couple of options for the subscription:

```cpp
subscribe(MessageType::Message, observer) // (1)!
subscribeAll(observer) // (2)!
subscribeAll({observer1, observer2}) // (3)!
unSubscribe(...) // (4)!
```

1. Subscribe to one specific message.
2. Subscribes to all the messages in `enum`.
3. Multiple observers can subscribe at once.
4. Unsubscribe from specific messages or all messages.

To see all available options for ``NonLinearSolverState``,
we refer to the file ``ikarus/solver/nonlinearsolver/solverstate.hh``.

To see all available options for ``ControlState``,
we refer to the file ``ikarus/controlroutines/controlstate.hh``.
