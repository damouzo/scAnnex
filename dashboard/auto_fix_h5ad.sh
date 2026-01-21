#!/bin/bash
#
# Auto-fix H5AD compatibility issues before launching dashboard
#

FIX_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
H5AD_PATH="$1"

if [[ -z "$H5AD_PATH" ]]; then
    echo "Usage: $0 <path/to/file.h5ad>"
    exit 1
fi

if [[ ! -f "$H5AD_PATH" ]]; then
    echo "Error: File not found: $H5AD_PATH"
    exit 1
fi

# Check if file has compatibility issues
echo "Checking H5AD compatibility..."

HAS_ISSUE=$(conda run -n scannex-dashboard python -c "
import h5py
import sys

try:
    with h5py.File('$H5AD_PATH', 'r') as f:
        if 'uns/log1p/base' in f:
            encoding_type = f['uns/log1p/base'].attrs.get('encoding-type', '')
            if encoding_type == b'null' or encoding_type == 'null':
                print('needs_fix')
                sys.exit(0)
        print('ok')
except Exception:
    print('ok')
" 2>/dev/null)

if [[ "$HAS_ISSUE" == "needs_fix" ]]; then
    echo "  ⚠ Compatibility issue detected - applying automatic fix..."
    conda run -n scannex-dashboard python "$FIX_SCRIPT_DIR/fix_h5ad_compatibility.py" "$H5AD_PATH" 2>&1 | grep -E "(Fixing|Creating|Found|Fixed|Backup)"
    echo "  ✓ Fix applied"
else
    echo "  ✓ File is compatible"
fi
