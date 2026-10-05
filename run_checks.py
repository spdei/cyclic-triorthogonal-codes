#!/usr/bin/env python3
"""Run every self-check (no pytest needed; `pytest tests/` also works)."""
import importlib, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tests'))
sys.path.insert(0, HERE)
modules = sorted(f[:-3] for f in os.listdir(os.path.join(HERE, 'tests')) if f.startswith('test_') and f.endswith('.py'))
if len(sys.argv) > 1:
    modules = [m for m in modules if any(a in m for a in sys.argv[1:])]
fails = 0
for mod in modules:
    M = importlib.import_module(mod)
    for name in sorted(n for n in dir(M) if n.startswith('test_')):
        t = time.time()
        try:
            getattr(M, name)()
            print(f"PASS  {mod}.{name}  ({time.time() - t:.1f}s)", flush=True)
        except Exception:
            fails += 1
            print(f"FAIL  {mod}.{name}", flush=True)
            traceback.print_exc()
print("\nALL CHECKS PASSED" if fails == 0 else f"\n{fails} CHECK(S) FAILED")
sys.exit(1 if fails else 0)
