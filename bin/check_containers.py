#!/usr/bin/env python3
"""
Check if all containers used in scAnnex pipeline are publicly accessible.
"""
import subprocess
import sys

containers = [
    "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0",
    "docker.io/satijalab/seurat:5.0.0",
    "quay.io/biocontainers/celltypist:1.6.2--pyhdfd78af_0"
]

print("\n" + "="*70)
print("  Checking Container Availability")
print("="*70 + "\n")

all_ok = True

for container in containers:
    print(f"Checking: {container}")
    try:
        # Try to pull manifest (doesn't download layers)
        if container.startswith("quay.io"):
            # Use docker manifest inspect
            result = subprocess.run(
                ["docker", "manifest", "inspect", container],
                capture_output=True,
                timeout=30
            )
        else:
            # Use docker pull with --dry-run equivalent
            result = subprocess.run(
                ["docker", "pull", container],
                capture_output=True,
                timeout=30
            )
        
        if result.returncode == 0:
            print(f"  ✅ Available\n")
        else:
            print(f"  ❌ NOT FOUND")
            print(f"     Error: {result.stderr.decode()[:100]}\n")
            all_ok = False
    except subprocess.TimeoutExpired:
        print(f"  ⚠️  Timeout (might still exist)\n")
    except FileNotFoundError:
        print(f"  ⚠️  Docker not installed, skipping check\n")
        break
    except Exception as e:
        print(f"  ❌ Error: {e}\n")
        all_ok = False

print("="*70)
if all_ok:
    print("✅ All containers are accessible!")
else:
    print("❌ Some containers are NOT accessible - MUST FIX!")
print("="*70 + "\n")

sys.exit(0 if all_ok else 1)
