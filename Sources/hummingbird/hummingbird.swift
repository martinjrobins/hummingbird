
import Numerics

protocol ButchersTableau {
  associatedtype Scalar
  var stages: Int { get }
  var order: Int { get }
  var estimator_order: Int { get }
  var dense_order: Int { get }
  var c: [Scalar] { get }
  var a: [[Scalar]] { get }
  var b_hat: [Scalar] { get }
  var b_star: [Scalar] { get }
  var b: [Scalar] { get }
  var p: [[Scalar]] { get }
}

struct DormondPrice<Scalar: Real & BinaryFloatingPoint>: ButchersTableau {
  typealias Scalar = Scalar
  let stages = 7
  let order = 5
  let estimator_order = 4
  let dense_order = 5
  let c: [Scalar] = 
    [Scalar(0.0), Scalar(1.0/5.0), Scalar(3.0/10.0), Scalar(4.0/5.0), Scalar(8.0/9.0),
      Scalar(1.0), Scalar(1.0)]
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
      [Scalar(35.0/384.0), Scalar(0.0), Scalar(500.0/1113.0), Scalar(125.0/192.0),
        Scalar(-2187.0/6784.0), Scalar(11.0/84.0)],
    ]
  let b_hat: [Scalar] = 
    [Scalar(35.0/384.0), Scalar(0.0), Scalar(500.0/1113.0), Scalar(125.0/192.0),
      Scalar(-2187.0/6784.0), Scalar(11.0/84.0)]
  let b: [Scalar] =
    [Scalar(5179.0/57600.0), Scalar(0.0), Scalar(7571.0/16695.0), Scalar(393.0/640.0),
      Scalar(-92097.0/339200.0), Scalar(187.0/2100.0), Scalar(1.0/40.0)]
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
}


func explicit_runge_kutta<Vector: HumVector, Tableau: ButchersTableau>(tableau: Tableau, 
                          ys: inout UnsafeMutableBufferPointer<Vector>, 
                          ts: Array<Vector.Scalar>, y0: Vector, 
                          dydx: @escaping (Vector, Vector.Scalar) -> Vector, 
                          tol: Vector.Scalar) -> Int where Vector.Scalar == Tableau.Scalar, Vector.Scalar: BinaryFloatingPoint {
                                                                            
  typealias Scalar = Tableau.Scalar
  let stages = tableau.stages
  let dense_order = tableau.dense_order
  let order = tableau.order
  let a = tableau.a
  let p = tableau.p
  let c = tableau.c
  let b = tableau.b
  let b_hat = tableau.b_hat
  assert(c.last! == 1.0, "last c value must be 1.0")

  var y_hat_n = y0
  ys[0] = y0
  var it = 1
  var k: [Vector] = Array<Vector>(repeating: Vector(repeating: 0), count: stages)

  let N = ts.count 
  if (N == 0) {
    return 0
  }
  var t_n = ts[0]
  var dydx_n = dydx(y0, t_n)
  var h_n = ts[N - 1] - t_n

  while (t_n < ts[N - 1]) { 
    var step_rejected = true 
    while (step_rejected) {
      // used to store the last dydx eval so can be re-used at the next step
      var dydx_stored = dydx_n

      // calculate stages
      k[0] = h_n * dydx_n
      var sum_ak: Vector     
      for i in 1..<stages {
        sum_ak = Vector(repeating: 0)     
        for j in 0..<i {
          sum_ak += a[i][j] * k[j]
        }
        // Note: if step is successful, this value at i=stages-1 
        // will be reused at the start of next step
        dydx_stored = dydx(y_hat_n + sum_ak, t_n + h_n * c[i])
        k[i] = h_n * dydx_stored
      }

      // calculate final value and error
      var error = Vector(repeating: 0)
      var y_hat_np1 = y_hat_n
      for i in 0..<stages {
        y_hat_np1 += b_hat[i] * k[i]
        error += (b_hat[i] - b[i]) * k[i]
      }

      let E_hp1 = error.inf_norm()
      if (E_hp1 < tol) {
        // if moved over any requested times then interpolate their values
        let t_np1 = t_n + h_n
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
        dydx_n = dydx_stored
        y_hat_n = y_hat_np1
        t_n = t_np1
      }
      
      // adapt step size
      h_n *= 0.9 * Scalar.pow(tol / E_hp1, 1.0/(Scalar(order) + 1.0))
    }
  }
  assert(it == ts.count)
  return it
}

extension HumVector where Scalar: BinaryFloatingPoint {

static func integrate(over ts: Array<Self.Scalar>, y0: Self, tol: Self.Scalar,
               dydx: @escaping (Self, Self.Scalar) -> Self) -> Array<Self> {
  let n = ts.count
  return Array<Self>(unsafeUninitializedCapacity: n) { buffer, initializedCount in
    initializedCount = explicit_runge_kutta(tableau: DormondPrice<Scalar>(), ys: &buffer, ts:
                                            ts, y0: y0, dydx: dydx, tol: tol)
  }
}

}

