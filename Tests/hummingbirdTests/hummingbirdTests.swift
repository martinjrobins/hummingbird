import XCTest
@testable import hummingbird

final class hummingbirdTests: XCTestCase {
    func testExample() {
        // This is an example of a functional test case.
        // Use XCTAssert and related functions to verify your tests produce the correct
        // results.
        XCTAssertEqual(hummingbird().text, "Hello, World!")
    }

    static var allTests = [
        ("testExample", testExample),
    ]
}
