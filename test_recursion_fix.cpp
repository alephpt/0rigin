// Simple test to verify operator splitting recursion fix
#include "src/SMFT.h"
#include "src/SMFTEngine.h"
#include <iostream>

int main() {
    std::cout << "=== Operator Splitting Recursion Fix Test ===" << std::endl;

    // Initialize Nova
    Nova* nova = new Nova();
    if (!nova->Init()) {
        std::cerr << "Failed to initialize Nova" << std::endl;
        return 1;
    }
    std::cout << "✓ Nova initialized" << std::endl;

    // Create SMFTEngine
    SMFTEngine engine(nova);
    std::cout << "✓ SMFTEngine created" << std::endl;

    // Initialize 64x64 grid
    if (!engine.initialize(64, 64, 2.5f, 0.0f)) {
        std::cerr << "Failed to initialize engine" << std::endl;
        return 1;
    }
    std::cout << "✓ Engine initialized (64x64)" << std::endl;

    // Initialize Dirac field
    engine.initializeDiracField(32.0f, 32.0f, 3.0f, 1.0f);
    std::cout << "✓ Dirac field initialized" << std::endl;

    // Test with N=2 (minimal operator splitting)
    std::cout << "\n--- Testing N=2 ---" << std::endl;
    try {
        engine.stepWithDirac(0.01f, 0.5f, 2, 2.0f, 0.05f);
        std::cout << "✓ N=2 test completed (no crash)" << std::endl;
    } catch (...) {
        std::cerr << "✗ N=2 test crashed!" << std::endl;
        return 1;
    }

    // Test with N=10
    std::cout << "\n--- Testing N=10 ---" << std::endl;
    try {
        engine.stepWithDirac(0.01f, 0.5f, 10, 2.0f, 0.05f);
        std::cout << "✓ N=10 test completed (no crash)" << std::endl;
    } catch (...) {
        std::cerr << "✗ N=10 test crashed!" << std::endl;
        return 1;
    }

    // Test with N=100
    std::cout << "\n--- Testing N=100 ---" << std::endl;
    try {
        engine.stepWithDirac(0.01f, 0.5f, 100, 2.0f, 0.05f);
        std::cout << "✓ N=100 test completed (no crash)" << std::endl;
    } catch (...) {
        std::cerr << "✗ N=100 test crashed!" << std::endl;
        return 1;
    }

    std::cout << "\n=== ALL TESTS PASSED ===" << std::endl;
    std::cout << "Recursion bug is FIXED!" << std::endl;

    delete nova;
    return 0;
}
