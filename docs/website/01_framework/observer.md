# Observer and Observable

To write output messages when desired by the user, the observer pattern is implemented in Ikarus.
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
  FINISHED_SUCESSFULLY,
  END
};
```

## IObservable

A class can be observable. The class then sends notifications when events are happening. To become observable, a class
must inherit from `IObservable<MessageType>`, for example,

```cpp
class NewtonRaphson : public IObservable<NonLinearSolverMessages> {...};
```

The function `this->notify(MessageType::Message)` is called at the appropriate position in the code to send a
notification. This could be, for example,

```cpp
this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
```

## IObserver

A class can be an observer. The class is then notified when events are happening and can perform actions. A very simple
example is shown below. To become an observer, the class must inherit from ``IObserver<MessageType>``, where ``MessageType``
is the `enum` of messages to use (see above).

```cpp
class OurFirstObserver : public IObserver<NonLinearSolverMessages> {
public:
  void updateImpl(NonLinearSolverMessages message) override {
    if (message == NonLinearSolverMessages::ITERATION_STARTED) std::cout << "Iteration started.\n";
  }
};
```

The observer has to implement the function ``void updateImpl(MessageType message)``. In this function, all actions can
be implemented that should be performed when the corresponding message is received.

To connect observer and observable, one has to call ``observalbe.subscribe(MessageType::Message,observer)``. Example:

```cpp
Ikarus::NewtonRaphson nr(...);
auto ourSimpleObserver = std::make_shared<OurFirstObserver>();
nr.subscribe(NonLinearSolverMessages::ITERATION_STARTED, ourSimpleObserver);
};
```

## Subscriptions

There are a couple of options for the subscription:

```cpp
subscribe(MessageType::Message,observer) // (1)!
subscribeAll(observer) // (2)!
subscribeAll({observer1,observer2}) // (3)!
unSubscribe(...) // (4)!
```

1. Subscribe to one specific message.
2. Subscribes to all the messages in `enum`.
3. Multiple observers can subscribe at once.
4. Unsubscribe from specific messages or all messages.

To send a message together with data, the sender (observable) calls

```cpp
this->notify(MessageType::Message, data);
```

and the receiver (observer) has to implement

```cpp
void updateImpl(MessageType message, data) override {...}
```

To see all available options for ``data``, we refer to the file ``observer.hh``.
