//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include "ikarus/utils/Observer/controlLogger.h"

class Control : public IObservable<ControlMessages> {
public:
  void solve() {}
};

TEST(Observer, ControlObserver) {
  auto controlObserver = std::make_shared<ControlLogger>();
  Control control;
  control.subscribeAll(controlObserver);

  control.notify(ControlMessages::CONTROL_STARTED);
  control.notify(ControlMessages::STEP_ENDED);
  control.notify(ControlMessages::SOLUTION_CHANGED);
}
