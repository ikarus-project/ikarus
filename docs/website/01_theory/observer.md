<!--
SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de

SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Oberserver and Observable

To solve situations like "I want to write an output when the loadstep is completed", the observer pattern is implemented in Ikarus.
To understand it and use it in your implementation, you need to understand three things: ``IObservable``, ``IObserver`` and Messages.

## Messages
A message class is a list of possible events that can happen and might be of interest. The messages which are used for 
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
A class can be observable. The class then sends notifications when events are happening. To become observable, a class has
to inherit from ``IObservable<MessageType>``, e.g.
```cpp
class NewtonRaphson : public IObservable<NonLinearSolverMessages> {...};
```
To send a notification, the function ``this->notify(MessageType::Message)`` is called at the corresponding position in the code.
This could be for example
```cpp
this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
```


## IObserver
A class can be an observer. The class is then notified when events are happening and can perform actions. A very simple
example is shown below. To become an observer, the class has to inherit from ``IObserver<MessageType>``, where ``MessageType``
is the enum of messages that should be used (see above). 
```cpp
class OurFirstObserver : public IObserver<NonLinearSolverMessages> {
public:
  void updateImpl(NonLinearSolverMessages message) override {
    if (message == NonLinearSolverMessages::ITERATION_STARTED) std::cout << "Yeah, the iteration started. Let's go!\n";
  }
};
```
The observer has to implement hat function ``void updateImpl(MessageType message)``. In this function, all actions can
be implemented that should be performed when the corresponding message is received.

To connect observer and observable, one has to call ``observalbe.subscribe(MessageType::Message,observer)``. Example:
```cpp
Ikarus::NewtonRaphson nr(...);
auto ourSimpleObserver = std::make_shared<OurFirstObserver>();
nr.subscribe(NonLinearSolverMessages::ITERATION_STARTED, ourSimpleObserver);
};
```

## Subscription options and sending data
There are a couple of options for the subscription:
```cpp
subscribe(MessageType::Message,observer) // (1)
subscribeAll(observer) // (2)
subscribeAll({observer1,observer2}) // (3)
unSubscribe(...) // (4)
```

1. Subscribe to one specific message.
2. Subscribe to all messages in the enum.
3. Multiple observers can subscribe at once.
4. Unsubscribe from specific messages or all messages

To send a message together with data, the sender (observable) calls
```cpp
this->notify(MessageType::Message,data);
```
and the receiver (observer) has to implement
```cpp
void updateImpl(MessageType message, data) override {
```
To see all available options for ``data``, we refer to the file ``observer.hh``.