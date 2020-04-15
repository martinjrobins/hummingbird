
func integrate<Vector: HumVector>(over ts: Array<Vector.Scalar>, y0: Vector, 
                              dydx: (Vector, Vector.Scalar) -> Vector) -> Array<Vector> {
  let dt = ts[1] - ts[0]
  let n = ts.count
  return Array<Vector>(unsafeUninitializedCapacity: n) { buffer, initializedCount in
    buffer[0] = y0
    for i in 1..<n {
        buffer[i] = buffer[i - 1] + dt * dydx(buffer[i - 1], ts[i - 1])
    }
    initializedCount = n
  }
}

