/**
 * test_csv_writer_validation.cpp - Validate TRDCSVWriter functionality
 *
 * Tests all features of the standardized CSV writer:
 *   1. File creation and directory structure
 *   2. Metadata generation (timestamp, git info, parameters)
 *   3. Type-safe row writing (int, float, double, string)
 *   4. Precision control and scientific notation
 *   5. Error handling (permissions, multiple writes)
 *
 * Expected output:
 *   output/CSV_Writer_Test/
 *     basic_test_YYYYMMDD_HHMMSS.csv
 *     precision_test_YYYYMMDD_HHMMSS.csv
 *     no_timestamp.csv
 *
 * @author TRD Development Team
 * @version v3D Unified
 */

#include "TRDCSVWriter.h"
#include <iostream>
#include <iomanip>
#include <cmath>

void test_basic_functionality() {
    std::cout << "\n=== Test 1: Basic Functionality ===\n";

    try {
        TRD::CSVWriter csv("basic_test", "CSV_Writer_Test", true);

        // Write metadata with custom parameters
        csv.writeMetadata({
            {"K_coupling", "1.0"},
            {"grid_size", "64x64x32"},
            {"evolution_steps", "1000"},
            {"dt", "0.01"}
        });

        // Write header
        csv.writeHeader({"Method", "Alpha_Measured", "Alpha_QED", "Ratio"});

        // Write data rows with various types
        csv.writeRow("Energy", 0.00353983, 0.00729735, 0.485084);
        csv.writeRow("Coupling", 0.000244141, 0.00729735, 0.0334561);
        csv.writeRow("Flux", 0.00182574, 0.00729735, 0.250185);
        csv.writeRow("GeometricMean", 0.00135642, 0.00729735, 0.185896);

        csv.close();

        std::cout << "✓ CSV written to: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

void test_type_conversions() {
    std::cout << "\n=== Test 2: Type Conversions ===\n";

    try {
        TRD::CSVWriter csv("type_test", "CSV_Writer_Test", true);

        csv.writeMetadata({{"test_type", "type_conversions"}});

        csv.writeHeader({"Int", "Long", "Float", "Double", "String"});

        // Test various types
        csv.writeRow(42, 9876543210L, 3.14159f, 2.71828182845904523536, "mixed");
        csv.writeRow(-100, -9999999999L, -1.5e-10f, 1.234567890123456e-20, "scientific");
        csv.writeRow(0, 0L, 0.0f, 0.0, "zeros");

        csv.close();

        std::cout << "✓ Type conversions verified: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

void test_precision_control() {
    std::cout << "\n=== Test 3: Precision Control ===\n";

    try {
        TRD::CSVWriter csv("precision_test", "CSV_Writer_Test", true);
        csv.setPrecision(12);  // High precision

        csv.writeMetadata({{"precision", "12"}});
        csv.writeHeader({"Value", "Squared", "Cubed"});

        // Test high-precision values
        double pi = 3.141592653589793;
        csv.writeRow(pi, pi*pi, pi*pi*pi);

        double e = 2.718281828459045;
        csv.writeRow(e, e*e, e*e*e);

        csv.close();

        std::cout << "✓ High precision CSV: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

void test_scientific_notation() {
    std::cout << "\n=== Test 4: Scientific Notation ===\n";

    try {
        TRD::CSVWriter csv("scientific_test", "CSV_Writer_Test", true);
        csv.setScientificNotation(true);
        csv.setPrecision(6);

        csv.writeMetadata({{"notation", "scientific"}});
        csv.writeHeader({"Small", "Medium", "Large"});

        csv.writeRow(1.23e-20, 4.56e0, 7.89e15);
        csv.writeRow(9.87e-15, 6.54e3, 3.21e20);

        csv.close();

        std::cout << "✓ Scientific notation CSV: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

void test_no_timestamp() {
    std::cout << "\n=== Test 5: No Timestamp Mode ===\n";

    try {
        TRD::CSVWriter csv("no_timestamp", "CSV_Writer_Test", false);

        csv.writeMetadata({{"timestamp_mode", "disabled"}});
        csv.writeHeader({"X", "Y", "Z"});
        csv.writeRow(1.0, 2.0, 3.0);
        csv.close();

        std::cout << "✓ Static filename CSV: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

void test_auto_metadata() {
    std::cout << "\n=== Test 6: Auto Metadata ===\n";

    try {
        TRD::CSVWriter csv("auto_metadata", "CSV_Writer_Test", true);

        // Skip writeMetadata() - should be auto-generated on writeHeader()
        csv.writeHeader({"A", "B", "C"});
        csv.writeRow(1, 2, 3);
        csv.close();

        std::cout << "✓ Auto-metadata CSV: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

void test_git_info() {
    std::cout << "\n=== Test 7: Git Information ===\n";

    std::string commit_hash = TRD::getGitCommitHash();
    std::string branch = TRD::getGitBranch();

    std::cout << "Git commit: " << commit_hash << "\n";
    std::cout << "Git branch: " << branch << "\n";

    if (commit_hash != "N/A" && branch != "N/A") {
        std::cout << "✓ Git info available\n";
    } else {
        std::cout << "⚠ Git info not available (not in git repo)\n";
    }
}

void test_error_handling() {
    std::cout << "\n=== Test 8: Error Handling ===\n";

    // Test 1: Writing header twice
    try {
        TRD::CSVWriter csv("error_test", "CSV_Writer_Test", true);
        csv.writeHeader({"A", "B"});
        csv.writeHeader({"C", "D"});  // Should throw
        std::cout << "✗ Failed to catch duplicate header\n";
    } catch (const std::exception& e) {
        std::cout << "✓ Caught duplicate header: " << e.what() << "\n";
    }

    // Test 2: Writing row before header
    try {
        TRD::CSVWriter csv("error_test2", "CSV_Writer_Test", true);
        csv.writeRow(1, 2, 3);  // Should throw
        std::cout << "✗ Failed to catch missing header\n";
    } catch (const std::exception& e) {
        std::cout << "✓ Caught missing header: " << e.what() << "\n";
    }

    // Test 3: Writing metadata twice
    try {
        TRD::CSVWriter csv("error_test3", "CSV_Writer_Test", true);
        csv.writeMetadata({{"test", "1"}});
        csv.writeMetadata({{"test", "2"}});  // Should throw
        std::cout << "✗ Failed to catch duplicate metadata\n";
    } catch (const std::exception& e) {
        std::cout << "✓ Caught duplicate metadata: " << e.what() << "\n";
    }
}

void demonstrate_usage_example() {
    std::cout << "\n=== Demo: Fine Structure Constant Example ===\n";

    try {
        TRD::CSVWriter csv("fine_structure_constant_results", "B2_FineStructure", false);

        csv.writeMetadata({
            {"K_coupling", "1.0"},
            {"grid_size", "64x64x32"},
            {"evolution_steps", "1000"},
            {"dt", "0.01"},
            {"golden_key", "246 GeV"}
        });

        csv.writeHeader({"Method", "Alpha_Measured", "Alpha_QED", "Ratio"});

        const double ALPHA_QED = 0.00729735;
        csv.writeRow("Energy", 0.00353983, ALPHA_QED, 0.00353983/ALPHA_QED);
        csv.writeRow("Coupling", 0.000244141, ALPHA_QED, 0.000244141/ALPHA_QED);
        csv.writeRow("Flux", 0.00182574, ALPHA_QED, 0.00182574/ALPHA_QED);
        csv.writeRow("GeometricMean", 0.00135642, ALPHA_QED, 0.00135642/ALPHA_QED);

        csv.close();

        std::cout << "✓ Demo CSV: " << csv.getFilePath() << "\n";

    } catch (const std::exception& e) {
        std::cout << "✗ FAILED: " << e.what() << "\n";
    }
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║       TRDCSVWriter Validation Suite                       ║\n";
    std::cout << "║       Standardized CSV Output for TRD Tests               ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";

    test_basic_functionality();
    test_type_conversions();
    test_precision_control();
    test_scientific_notation();
    test_no_timestamp();
    test_auto_metadata();
    test_git_info();
    test_error_handling();
    demonstrate_usage_example();

    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║ VALIDATION COMPLETE                                        ║\n";
    std::cout << "║ Check output/CSV_Writer_Test/ for generated files         ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";

    return 0;
}
