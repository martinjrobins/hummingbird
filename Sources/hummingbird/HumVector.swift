import Numerics

public protocol HumVector: AdditiveArithmetic {
  associatedtype Scalar: Real, BinaryFloatingPoint

  var scalarCount: Int { get }
  
  subscript(index: Int) -> Self.Scalar { get set }

  init(repeating x: Scalar)

  static func * (lhs: Self.Scalar, rhs: Self) -> Self

  func inf_norm() -> Scalar 
}

// default implementation for inf_norm
extension HumVector {

  public func inf_norm() -> Scalar { 
    var max: Scalar = abs(self[0])
    for i in 1..<self.scalarCount {
      let abs_data = abs(self[i])
      if (abs_data > max) {
        max = abs_data
      }
    }
    return max
  }

}

// make a Double conform to a HumVector of size 1
extension Double : HumVector {
  public typealias Scalar = Double

  public subscript(index: Int) -> Self.Scalar {
    get {
      return self
    } 
    set {
      self = newValue
    } 
  }

  public var scalarCount: Int { return 1 }

  public init(repeating x: Scalar) { self = x }

  public func inf_norm() -> Scalar { return abs(self) } 
}

// make a Float conform to a HumVector of size 1
extension Float: HumVector {
  public typealias Scalar = Float

  public subscript(index: Int) -> Self.Scalar {
    get {
      return self
    } 
    set {
      self = newValue
    } 
  }

  public var scalarCount: Int { return 1 }

  public init(repeating x: Scalar) { self = x }

  public func inf_norm() -> Scalar { return abs(self) } 
}

// make all the SIMD types conform to a HumVector
extension SIMD2: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD2: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD3: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD3: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD4: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD4: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD8: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD8: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD16: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD16: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD32: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD32: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD64: HumVector where Scalar: Real & BinaryFloatingPoint {
}

extension SIMD64: AdditiveArithmetic where Scalar: Real & BinaryFloatingPoint {
}
