struct DormondPriceStepper<Vector: HumVector> {
  var y_hat_n: Vector
  var t_n: Vector.Scalar
  var it: Int = 0
  var ys: Array<Vector>
  let ts: Array<Vector.Scalar>
  let dydx: (Vector, Vector.Scalar) -> Vector
  var k: Array<Vector>
  let tol: Vector.Scalar
  let s = 7
  let p = 5
  let q = 4
  let a: [[Vector.Scalar]] = [
    [0.0,             0.0,            0.0,            0.0,          0.0,             0.0      ], 
    [1.0/5.0,         0.0,            0.0,            0.0,          0.0,             0.0      ], 
    [3.0/40.0,        9.0/40.0,       0.0,            0.0,          0.0,             0.0      ], 
    [44.0/45.0,      -56.0/15.0,      32.0/9.0,       0.0,          0.0,             0.0      ], 
    [19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0, -212.0/729.0, 0.0,             0.0      ], 
    [-9017.0/3168.0, -355.0/33.0,     46732.0/5247.0, 49.0/176.0,   -5103.0/18656.0, 0.0      ], 
    [35.0/386.0,      0.0,            500.0/1113.0,   125.0/192.0,  -2187.0/6784.0,  11.0/84.0], 
  ]
  let b_hat: [Vector.Scalar] = 
    [35.0/384.0, 0.0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0, 0.0]
  //let b: [Vector.Scalar] = 
  //  [5179.0/57600.0, 0.0, 7571.0/16695.0, 393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0]
  let b: [Vector.Scalar] = 
    [1951.0/21600.0, 0.0, 22642.0/50085.0, 451.0/720.0, -12231.0/42400.0, 649.0/6300.0, 1.0/60.0]
  let b_star: [Vector.Scalar] = 
    [6025192743.0/30085553152.0, 0.0, 51252292925.0/65400821598.0,
      -2691868925.0/45128329728.0, 187940372067.0/1594534317056.0,
      -1776094331.0/19743644256.0, 11237099.0/235043384.0]


  init(ys: inout Array<Vector>, ts: Array<Vector.Scalar>, y0: Vector, 
       dydx: @escaping (Vector, Vector.Scalar) -> Vector, tol: Vector.Scalar) {
    self.y_hat_n = y0
    self.ts = ts
    self.dydx = dydx
    self.t_n = ts[0]
    self.ys[0] = y0
    self.it = 1
    self.ys = ys;
    self.k = [Vector](repeating: Vector(repeating: 0), count: s)
    self.tol = tol
  }

  func step() {
    var step_rejected = true 
    guard let h_n = ts.last - t_n else { return }
    while (step_rejected) {
      k[0] = h_n * dydx(y_hat_n)
      for i in 1..<s {
        var sum_ak = Vector(repeating: 0)     
        for j in 0..<i {
          sum_ak += a[i][j] * k[j]
        }
        k[i] = h_n * dydx(y_hat_n + sum_ak)
      }
      var y_np1 = y_hat_n
      var y_hat_np1 = y_hat_n
      for i in 0..<s {
        y_np1 += b[i] * k[i]
        y_hat_np1 += b_hat[i] * k[i]
      }
      var E_hp1 = (y_hat_np1 - y_np1).inf_norm()
      if (E_hp1 < tol) {
        step_rejected = false
        let t_np1 = t_n + h_n
        if (t_n > ts[it]) {
          let t_nph = t_n + 0.5 * h_n
          var y_nph = y_hat_n
          for i in 0..<s {
            y_nph += 0.5 * b_star[i] * k[i]
          }
          let interpolator = QuarticInterpolation(y_hat_n, y_nph, y_hat_np1, 
                                            t_n, t_nph, t_np1) 
          do {
            ys[it] = interpolator(ts[it])
            it += 1
          } while (t_n > ts[it])
        }
      }
      h_n = 0.9 * h_n pow(self.tol / E_hp1, 1.0/(p + 1))
    }
  }

}


func dormond_price_step<Vector: HumVector>(y0: Vector, from t0: Vector.Scalar, to t1:
                                           Vector.Scalar, 
                              dydx: (Vector, Vector.Scalar) -> Vector) -> Array<Vector> {
  let dt = ts[1] - ts[0]
  let n = ts.count
  return Array<Vector>(unsafeUninitializedCapacity: n) { buffer, initializedCount in
    buffer[0] = y0
    var stepper = DormondPriceStepper(ts, y0, dydx); 
    while (stepper.t0 < t0.last) {
      stepper.step();
    }
    initializedCount = n
  }
}

func dormond_price<Vector: HumVector>(over ts: Array<Vector.Scalar>, y0: Vector, 
                              dydx: (Vector, Vector.Scalar) -> Vector) -> Array<Vector> {
  let dt = ts[1] - ts[0]
  let n = ts.count
  return Array<Vector>(unsafeUninitializedCapacity: n) { buffer, initializedCount in
    buffer[0] = y0
    for i in 1..<n {
      dormond_price_step(y0: buffer[i - 1], from: ts[i - 1], to: ts[i], dydx: dydx)
    }
    initializedCount = n
  }
}

