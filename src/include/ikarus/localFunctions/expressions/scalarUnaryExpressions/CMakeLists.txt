# install headers
install(
  FILES scalarUnaryExpression.hh scalarSqrtExpr.hh scalarPowExpr.hh scalarLogExpr.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/localFunctions/expressions/scalarUnaryExpressions
)
