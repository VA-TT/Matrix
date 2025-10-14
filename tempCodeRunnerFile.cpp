  Vector<Matrix<double, d, d * nNodes>> C(nNodes);

  // Build connectivity matrices C[n] (d x (d * nNodes)) where an identity
  // block of size dxd is placed at the columns corresponding to node n.
  // Note: d and nNodes are constexpr so d * nNodes is a compile-time
  // constant and the Matrix template can be instantiated.
  const auto Id = Matrix<double, d, d>::identity();
  for (Index n{0}; n < nNodes; ++n) {
    Matrix<double, d, d * nNodes> Ci{}; // initialized to zero
    for (Index i = 0; i < d; ++i) {
      for (Index j = 0; j < d; ++j) {
        Ci(i, n * d + j) = Id(i, j);
      }
    }
    C[n] = Ci;
  }
  std::cout << "Checkpoint: constructed C matrices" << std::endl;
  std::cout << C;