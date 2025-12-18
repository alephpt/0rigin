# Security Fix: Command Injection Vulnerabilities Resolved

## Summary
Fixed 4 command injection vulnerabilities by replacing `system()` calls with `std::filesystem::create_directories()`.

## Files Modified

### 1. `test/test_stochastic_cpu.cpp`
- **Line 68**: Replaced `system("mkdir -p output/stochastic_cpu")` with `std::filesystem::create_directories("/home/persist/neotec/0rigin/output")`
- Added `#include <filesystem>` header
- Removed `#include <sys/stat.h>` header (no longer needed)

### 2. `test/test_smft_headless.cpp`
- **Line 343**: Replaced `system("mkdir -p build/output/headless")` with `std::filesystem::create_directories("output")`
- Added `#include <filesystem>` header

### 3. `test/test_noise_sweep_corrected.cpp`
- **Line 95**: Replaced `system("mkdir -p output/noise_sweep")` with `std::filesystem::create_directories("/home/persist/neotec/0rigin/output/corrected_sweep")`
- Added `#include <filesystem>` header

### 4. `test/test_output_structure.cpp`
- **Line 49**: Replaced `system(create_cmd.c_str())` with `std::filesystem::create_directories(dirname)`
- Added `#include <filesystem>` header
- Removed `#include <sys/stat.h>` and `#include <sys/types.h>` headers (no longer needed)
- Refactored from using char arrays to std::string for directory and file paths

## Security Impact
- **Before**: `system()` calls could potentially execute arbitrary shell commands if user input was incorporated into paths
- **After**: `std::filesystem::create_directories()` safely creates directories without invoking a shell, eliminating command injection risk

## Technical Details
- Uses C++17 `std::filesystem` library for safe, cross-platform directory creation
- `create_directories()` creates parent directories as needed (equivalent to `mkdir -p`)
- Returns without error if directory already exists
- Throws `std::filesystem::filesystem_error` on failure (can be caught if needed)

## Compilation Requirements
- Requires C++17 or later (`-std=c++17` flag)
- May need to link with `-lstdc++fs` on some older compilers

## Verification
```bash
# Verify no system() calls remain
grep -n "system(" test/test_stochastic_cpu.cpp test/test_smft_headless.cpp \
                   test/test_noise_sweep_corrected.cpp test/test_output_structure.cpp
# Should return no results

# Verify filesystem usage
grep -n "std::filesystem::create_directories" test/*.cpp
# Should show 4 occurrences
```

## Result
✅ All 4 command injection vulnerabilities have been successfully eliminated.