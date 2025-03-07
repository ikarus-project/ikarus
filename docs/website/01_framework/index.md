# Framework

``` mermaid
classDiagram
  GridView <-- Grid
  GlobalBasis <-- GridView
  Assembler <-- GlobalBasis
  DifferentiableFunction <-- Assembler
  Assembler <-- FiniteElement
  FiniteElement <-- FERequirements
  FERequirements <|-- ResultRequirements
  FiniteElement <-- ResultRequirements
  FiniteElement <-- Local function
  Localfunction <-- Localbasis
  GlobalBasis <--> Localbasis
  Localfunction <-- Manifold
  NonlinearSolver <-- DifferentiableFunction
  NonlinearSolver <-- LinearSolver
  Controlroutine <-- NonlinearSolver
  VTKWriter <-- Controlroutine
  DirichletConditions .. Assembler
  DirichletConditions .. Controlroutine
  Controlroutine <|-- IObservable
  NonlinearSolver <|-- IObservable
  Observer ..> IObservable
  FERequirements <-- Affordances
  Affordances <-- ScalarAffordance
  Affordances <-- VectorAffordance
  Affordances <-- MatrixAffordance
  ResultRequirements <-- ResultType


  class ScalarAffordance{
  <<enumeration>>
      mechanicalPotentialEnergy
      microMagneticPotentialEnergy
      ...
  }

  class VectorAffordance{
        <<enumeration>>
      forces
      microMagneticForces
      ...
  }

  class MatrixAffordance{
        <<enumeration>>
      stiffness
      materialstiffness
      geometricstiffness
      mass
      stiffnessdiffBucklingVector
      microMagneticHessian
      ...
  }

  class ResultType{
      <<enumeration>>
      noType
      magnetization
      gradientNormOfMagnetization
      vectorPotential
      divergenceOfVectorPotential
      BField
      HField
      cauchyStress
      director
      ...
  }

  class Observer{
    +update()
  }

  class DirichletConditions{
  TBA
  }

  class IObservable{
    +subscribe()
    +subscribeAll()
    +unSubscribe()
    +unSubscribeAll()
    +notify()
  }

  class FERequirements{
    +hasAffordance()
    +getGlobalSolution()
    +getParameter()
  }

  class ResultRequirements{
    +isResultRequested()
    +getParameter()
  }

  class Assembler{
    +scalar()
    +matrix()
    +vector()
    +createFullVector()
  }
  class GridView{
    +elements(gridView)
    +vertices(gridView)
    +edges(gridView)
    +surfaces(gridView)
  }
  class Controlroutine{
    +run()
  }

  class NonlinearSolver{
    +setup()
    +solve()
    +differentiableFunction()
  }

  class GlobalBasis{
    +localView()
  }

    class Grid{
    +leafGridView()
  }

  class FiniteElement{
  +calculateScalar()
  +calculateVector()
  +calculateMatrix()
  +calculateAt()
  }

  class LinearSolver{
    +analyzePattern()
    +factorize()
    +compute()
    +solve()
  }
  class DifferentiableFunction{
    +operator()
    +derivative()
}

  class Localfunction{
    +calculateFunction()
    +calculateDerivative()
    +bind()
    +viewOverIntegrationPoints()
  }

  class Localbasis{
    +calculateFunction()
    +evaluateJacobian()
    +bind()
    +isBound()
    +viewOverIntegrationPoints()
  }

  class Manifold{
    +setValue()
    +operator+=()
    +getValue()
    +size()
    +size()
  }

click DifferentiableFunction href "../differentiablefunction/"
click LinearSolver href "../solvers/#linear-solver"
click NonlinearSolver href "../solvers/#non-linear-solver"
click FiniteElement href "../finiteElements/"
click GridView href "../grid/"
click Grid href "../grid/"
click Controlroutine href "../controlRoutines/"
click Assembler href "../assembler/"
click Localfunction href "../localFunctions/"
click Manifold href "../manifolds/"
click Localbasis href "../localBasis/"
click FERequirements href "../feRequirements/"
click ResultRequirements href "../feRequirements/#fe-result-requirements"
click Affordances href "../feRequirements/"
click ResultType href "../feRequirements/"
click IObservable href "../observer/#iobservable"
click Observer href "../observer/#iobserver"
click GlobalBasis href "../globalBasis/"

```
