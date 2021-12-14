//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <ikarus/utils/controlLogger.h>

class Control : public IObservable<ControlMessages>
{
public:
  void solve()
  {}

};

TEST(Observer, ControlObserver) {

  auto controlObserver = std::make_shared<ControlLogger>();
  Control control;
  control.subscribeAll(controlObserver);

control.notify(ControlMessages::LOADSTEP_ENDED);
control.notify(ControlMessages::ITERATION_ENDED);
control.notify(ControlMessages::RESIDUALNORM_UPDATED,5.0);

}


