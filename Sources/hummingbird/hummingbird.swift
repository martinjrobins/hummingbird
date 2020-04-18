
import Numerics

class DormondPriceStepper<Vector: HumVector> {
  typealias Scalar = Vector.Scalar
  var y_hat_n: Vector
  var t_n: Scalar
  var h_n: Scalar
  var dydx_n: Vector
  var it: Int = 0
  var ys: UnsafeMutableBufferPointer<Vector>
  let ts: Array<Scalar>
  let dydx: (Vector, Scalar) -> Vector
  var k: Array<Vector>
  let tol: Scalar
  let stages = 6
  let order = 5
  let estimator_order = 4
  let dense_order = 5
  let c: [Scalar] = 
    [Scalar(0.0), Scalar(1.0/5.0), Scalar(3.0/10.0), Scalar(4.0/5.0), Scalar(8.0/9.0),
      Scalar(1.0)]
  //let a: [[Scalar]] = 
  //  [
  //    [0.0,             0.0,            0.0,
  //      0.0,
  //      0.0], 
  //    [1.0/5.0,         0.0,            0.0,
  //      0.0,          0.0], 
  //    [3.0/40.0,        9.0/40.0,       0.0,
  //      0.0,          0.0], 
  //    [44.0/45.0,      -56.0/15.0,
  //      32.0/9.0,       0.0,          0.0], 
  //    [19372.0/6561.0, -25360.0/2187.0,
  //      64448.0/6561.0, -212.0/729.0, 0.0], 
  //    [-9017.0/3168.0, -355.0/33.0,
  //      46732.0/5247.0, 49.0/176.0,
  //      -5103.0/18656.0],
  //  ]
  let a: [[Scalar]] = 
    [
      [Scalar(0.0),             Scalar(0.0),            Scalar(0.0),
        Scalar(0.0),
        Scalar(0.0)], 
      [Scalar(1.0/5.0),         Scalar(0.0),            Scalar(0.0),
        Scalar(0.0),          Scalar(0.0)], 
      [Scalar(3.0/40.0),        Scalar(9.0/40.0),       Scalar(0.0),
        Scalar(0.0),          Scalar(0.0)], 
      [Scalar(44.0/45.0),      Scalar(-56.0/15.0),
        Scalar(32.0/9.0),       Scalar(0.0),          Scalar(0.0)], 
      [Scalar(19372.0/6561.0), Scalar(-25360.0/2187.0),
        Scalar(64448.0/6561.0), Scalar(-212.0/729.0), Scalar(0.0)], 
      [Scalar(-9017.0/3168.0), Scalar(-355.0/33.0),
        Scalar(46732.0/5247.0), Scalar(49.0/176.0),
        Scalar(-5103.0/18656.0)],
    ]
  let b_hat: [Scalar] = 
    [Scalar(35.0/384.0), Scalar(0.0), Scalar(500.0/1113.0), Scalar(125.0/192.0),
      Scalar(-2187.0/6784.0), Scalar(11.0/84.0)]
  let b: [Scalar] =
    [Scalar(5179.0/57600.0), Scalar(0.0), Scalar(7571.0/16695.0), Scalar(393.0/640.0),
      Scalar(-92097.0/339200.0), Scalar(187.0/2100.0), Scalar(1.0/40.0)]
  //let b: [Scalar] = 
  //  [1951.0/21600.0, 0.0, 22642.0/50085.0, 451.0/720.0, -12231.0/42400.0, 649.0/6300.0, 1.0/60.0]
  let b_star: [Scalar] =
    [Scalar(6025192743.0/30085553152.0), Scalar(0.0),
      Scalar(51252292925.0/65400821598.0),
      Scalar(-2691868925.0/45128329728.0), Scalar(187940372067.0/1594534317056.0),
      Scalar(-1776094331.0/19743644256.0), Scalar(11237099.0/235043384.0)]

  let p: [[Scalar]] = 
    [
      [Scalar(1.0), Scalar(-32272833064.0/11282082432.0),
        Scalar(34969693132.0/11282082432.0),
        Scalar(-13107642775.0/11282082432.0),
        Scalar(157015080.0/11282082432.0)],
      [Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0), Scalar(0.0)],
      [Scalar(0.0), Scalar(1323431896.0*100.0/32700410799.0),
        Scalar(-2074956840.0*100.0/32700410799.0),
        Scalar(914128567.0*100.0/32700410799.0),
        Scalar(-15701508.0*100.0/32700410799.0)],
      [Scalar(0.0), Scalar(-889289856.0*25.0/5641041216.0),
        Scalar(2460397220.0*25.0/5641041216.0),
        Scalar(-1518414297.0*25.0/5641041216.0),
        Scalar(94209048.0*25.0/5641041216.0)],
      [Scalar(0.0), Scalar(259006536.0*2187.0/199316789632.0),
        Scalar(-687873124.0*2187.0/199316789632.0),
        Scalar(451824525.0*2187.0/199316789632.0),
        Scalar(-52338360.0*2187.0/199316789632.0)], 
      [Scalar(0.0), Scalar(-361440756.0*11.0/2467955532.0),
        Scalar(946554244.0*11.0/2467955532.0),
        Scalar(-661884105.0*11.0/2467955532.0),
        Scalar(106151040.0*11.0/2467955532.0)],
      [Scalar(0.0), Scalar(44764047.0/29380423.0), Scalar(-127201567/29380423.0), 
        Scalar(90730570.0/29380423.0), Scalar(-8293050.0/29380423.0)],
    ]

   //P = np.array([
   //     [1, -8048581381/2820520608, 8663915743/2820520608,
   //      -12715105075/11282082432],
   //     [0, 0, 0, 0],
   //     [0, 131558114200/32700410799, -68118460800/10900136933,
   //      87487479700/32700410799],
   //     [0, -1754552775/470086768, 14199869525/1410260304,
   //      -10690763975/1880347072],
   //     [0, 127303824393/49829197408, -318862633887/49829197408,
   //      701980252875 / 199316789632],
   //     [0, -282668133/205662961, 2019193451/616988883, -1453857185/822651844],
   //     [0, 40617522/29380423, -110615467/29380423, 69997945/29380423]])


  init(ys: inout UnsafeMutableBufferPointer<Vector>, ts: Array<Scalar>, y0: Vector, 
       dydx: @escaping (Vector, Scalar) -> Vector, tol: Scalar) {
    self.y_hat_n = y0
    self.ts = ts
    self.dydx = dydx
    self.ys = ys;
    self.ys[0] = y0
    self.it = 1
    self.k = [Vector](repeating: Vector(repeating: 0), count: stages + 1)
    self.tol = tol

    let n = ts.count 
    if (n == 0) {
      self.t_n = Scalar(0.0)
      self.dydx_n = dydx(y0, self.t_n)
      self.h_n = Scalar(0.0)
      return
    }
    self.t_n = ts[0]
    self.dydx_n = dydx(y0, self.t_n)
    self.h_n = ts[n - 1] - self.t_n
  }

  func step() -> Scalar {
    var step_rejected = true 
    while (step_rejected) {
      k[0] = h_n * dydx_n
      for i in 1..<stages {
        var sum_ak = Vector(repeating: 0)     
        for j in 0..<i {
          sum_ak += a[i][j] * k[j]
        }
        k[i] = h_n * dydx(y_hat_n + sum_ak, t_n + h_n * c[i])
      }
      var sum_ak_final = Vector(repeating: 0)
      var y_np1 = y_hat_n
      // Note: b_hat is both the:
      //        - last row of a[i,j]
      //        - b values for the 5th order method
      for i in 0..<stages {
        y_np1 += b[i] * k[i]
        sum_ak_final += b_hat[i] * k[i]
      }
      // so can calculate the next position for the 5th order method using sum_ak_final
      let y_hat_np1 = sum_ak_final + y_hat_n  

      // and then eval the rhs using this position
      let t_np1 = t_n + h_n
      let dydx_np1 = dydx(y_hat_np1, t_np1)

      // then use the final k to calculate the next value for the 4th order
      // approximation
      k[stages] = h_n * dydx_np1
      y_np1 += b[stages] * k[stages]

      let E_hp1 = (y_hat_np1 - y_np1).inf_norm()
      if (E_hp1 < tol) {

        // if moved over any requested times then interpolate their values
        while (it < ts.count && t_np1 >= ts[it]) {
          let sigma = (ts[it] - t_n) / h_n
          var sigma_i = sigma
          var y_sigma = y_hat_n
          for i in 1...stages {
            sigma_i = sigma 
            var b_i = sigma_i * p[i][0]
            for j in 1..<dense_order {
              sigma_i *= sigma
              b_i += sigma_i * p[i][j]
            }
            y_sigma += b_i * k[i]
          }
          ys[it] = y_sigma
          it += 1
        }

        // move to next step
        step_rejected = false
        dydx_n = dydx_np1
        y_hat_n = y_hat_np1
        t_n = t_np1
      }
      
      // adapt step size
      h_n *= 0.9 * Scalar.pow(self.tol / E_hp1, 1.0/(Scalar(order) + 1.0))
    }
    return t_n
  }
}

func integrate<Vector: HumVector>(over ts: Array<Vector.Scalar>, y0: Vector, tol:
                                  Vector.Scalar,
                              dydx: @escaping (Vector, Vector.Scalar) -> Vector) -> Array<Vector> {
  let n = ts.count
  return Array<Vector>(unsafeUninitializedCapacity: n) { buffer, initializedCount in
    guard let ts_last = ts.last else { 
      initializedCount = 0
      return 
    }
    let stepper = DormondPriceStepper(ys: &buffer, ts: ts, y0: y0, dydx: dydx, tol: tol)
    while (stepper.step() < ts_last) {}
    initializedCount = n
  }
}

