"""Run all five supplied Python examples from any working directory."""
from pathlib import Path
import subprocess
import sys
root=Path(__file__).resolve().parent
for tier in range(1,6):
    print(f"\nRunning Tier {tier}",flush=True)
    subprocess.run([sys.executable,str(root/f"advanced/tier-{tier}/python/run_case.py")],cwd=root,check=True)
print("All example runs completed. Results are inside each tier's results/python folder.")
