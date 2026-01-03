#!/bin/bash
# Comprehensive SMFT → TRD Framework Rename Script
# Renames: Synchronization Mass Field Theory → Topological Resonance Dynamics
#          SMFT → TRD

set -e

echo "===== TRD Framework Rename Script ====="
echo "Renaming SMFT → TRD across entire codebase"
echo ""

# Function to perform safe replacement in a file
replace_in_file() {
    local file="$1"
    if [ -f "$file" ]; then
        # Create backup
        cp "$file" "$file.bak"

        # Perform replacements (order matters - specific to general)
        sed -i 's/Synchronization Mass Field Theory/Topological Resonance Dynamics/g' "$file"
        sed -i 's/SMFTEngine/TRDEngine/g' "$file"
        sed -i 's/SMFTCore/TRDCore/g' "$file"
        sed -i 's/SMFTBuffer/TRDBuffer/g' "$file"
        sed -i 's/SMFTCompute/TRDCompute/g' "$file"
        sed -i 's/SMFTDescriptor/TRDDescriptor/g' "$file"
        sed -i 's/SMFTPipeline/TRDPipeline/g' "$file"
        sed -i 's/SMFTTest/TRDTest/g' "$file"
        sed -i 's/MSFTCommon/TRDCommon/g' "$file"
        sed -i 's/MSFTEngine/TRDEngine/g' "$file"
        sed -i 's/MSFT/TRD/g' "$file"
        sed -i 's/SMFT/TRD/g' "$file"
        sed -i 's/smft/trd/g' "$file"

        # Check if file changed
        if ! diff -q "$file" "$file.bak" > /dev/null 2>&1; then
            echo "  ✓ Updated: $file"
            rm "$file.bak"
            return 0
        else
            # No changes, restore original
            mv "$file.bak" "$file"
            return 1
        fi
    fi
}

# Counter
updated=0

# 1. Update C++ source files
echo "[1/7] Updating C++ source files..."
for file in $(find . -type f \( -name "*.cpp" -o -name "*.h" \) ! -path "*/build/*" ! -path "*/.git/*" ! -path "*/lib/*"); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

# 2. Update YAML configs
echo "[2/7] Updating YAML configuration files..."
for file in $(find . -type f -name "*.yaml" ! -path "*/build/*" ! -path "*/.git/*"); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

# 3. Update Markdown docs
echo "[3/7] Updating Markdown documentation..."
for file in $(find . -type f -name "*.md" ! -path "*/build/*" ! -path "*/.git/*" ! -path "*/lib/*"); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

# 4. Update CMakeLists
echo "[4/7] Updating CMake files..."
for file in $(find . -type f -name "CMakeLists.txt" ! -path "*/build/*" ! -path "*/.git/*" ! -path "*/lib/*"); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

# 5. Update shader files
echo "[5/7] Updating shader files..."
for file in $(find . -type f -name "*.comp" ! -path "*/build/*" ! -path "*/.git/*"); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

# 6. Update scripts
echo "[6/7] Updating script files..."
for file in $(find ./scripts -type f -name "*.sh" 2>/dev/null); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

# 7. Update text files
echo "[7/7] Updating text files..."
for file in $(find . -type f -name "*.txt" ! -path "*/build/*" ! -path "*/.git/*" ! -path "*/lib/*"); do
    if replace_in_file "$file"; then
        ((updated++))
    fi
done

echo ""
echo "===== Summary ====="
echo "Files updated: $updated"
echo ""
echo "===== File Renames Needed (Manual) ====="
echo "The following files should be renamed:"
find . -type f \( -name "*SMFT*" -o -name "*smft*" -o -name "*MSFT*" \) ! -path "*/build/*" ! -path "*/.git/*" ! -path "*/lib/*" | head -20
echo ""
echo "Run this script again after renaming files to update references."
echo "===== Done ====="
