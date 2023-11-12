# Manifold elements

Usually optimization problems are defined in terms of a cost function, such as:

$$
 \min_{\boldsymbol{x} \in \mathcal{M}} f(\boldsymbol{x} )
$$
where \( f: \mathcal{M} \rightarrow \mathbb{R} \).

Usually \( \mathcal{M} \) is an Euclidean vector space \( \mathbb{R}^n \).

In a finite element context, for example, in 2D elasticity problems, we have a
two-dimensional displacement for each node. As a result, if there are \(n \) nodes, we optimize in \( {\mathbb{R}^{2n}} \).

The nodal degrees of freedom should be wrapped in `#!cpp Ikarus::RealTuple<double,2>` in this case.

Another case of optimization is on non-linear manifolds. These arise typically for Cosserat materials \( \mathcal{S}\mathcal{O}(3) \)[@sanderBeam],
Reissner-Mindlin shells \( \mathcal{S}^2 \)[@muller2022consistent] and micro-magnetics \( \mathcal{S}^{2} \) or incompressible materials.

## Interface

The general interface of the manifold elements is represented by the following concept.

```cpp
namespace Ikarus::Concepts {
  template <typename ManifoldType>
  concept Manifold = requires(ManifoldType var, typename ManifoldType::CorrectionType correction, std::ostream& s,
                              typename ManifoldType::CoordinateType value) {
    typename ManifoldType::ctype; // (1)!
    ManifoldType::valueSize; // (2)!
    ManifoldType::correctionSize; // (3)!
    typename ManifoldType::CoordinateType; // (4)!
    typename ManifoldType::CorrectionType; // (5)!
    { var.getValue() } -> std::convertible_to<typename ManifoldType::CoordinateType>; // (6)!
    { var.setValue(value) } -> std::same_as<void>; // (7)!
    { var+=correction };  // (8)!
    //...
  };
}
```

1. The type for the coordinate values, usually `#!cpp double`.
2. The number of values to store for the state of the element. E.g., the three-dimensional unit vector needs three entries to store its
   state.
3. The size of the correction for an element. `valueSize` and `correctionSize` are the same in Euclidean space. But, for example, the
   three-dimensional unit vector needs a two-dimensional correction (which lives in the tangent space).
4. The type to store the element coordinates is usually `#!cpp Eigen::Vector<double,ManifoldType::valueSize>`
5. The type to store the element correction type is usually `#!cpp Eigen::Vector<double,ManifoldType::correctionSize>`
6. Access the underlying coordinate vector of the manifold element.
7. Directly set the value. E.g., set `#!cpp Ikarus::UnitVector<double,3> a; a.setValue(Eigen::Vector3d::UnitZ());`
8. Update the element with a correction vector. E.g.,

     ```cpp
        Ikarus::RealTuple<double,3> a;
        a+= Eigen::Vector3d::UnitX();

        Ikarus::UnitVector<double,3> b;
        b+= Eigen::Vector2d::UnitX();
     ```

## Implementations

| Name                                 | Formal definition                                                                                                                                                     | Notes | Header          |
|:-------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------|-----------------|
| \(n\)-th dimensional Euclidean space | $$ \boldsymbol{x} \in \mathbb{R}^n  $$                                                                                                                                |       | `realTuple.hh`  |
| Unit sphere                          | $$ \boldsymbol{x} \in \mathcal{S}^{n-1}, \quad \mathcal{S}^{n-1} = \left\\{ \boldsymbol{x} \in \mathbb{R}^n :  \boldsymbol{x}\cdot  \boldsymbol{x}  = 1 \right\\}  $$ |       | `unitVector.hh` |

\bibliography
