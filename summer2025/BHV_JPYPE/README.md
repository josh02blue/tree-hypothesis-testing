## BHV distance via JPype (`gtp.jar`)

This module wraps the Java geodesic BHV distance implementation in `gtp.jar`
and exposes a simple Python function:

- `compute_bhv_distance(tree1, tree2, rooted=True) -> float`

### Requirements

- Python 3.x
- Java (JRE/JDK) on PATH
- `pip install jpype1`
- `gtp.jar` present; update `GTP_JAR` in `bhv_jpype.py` if needed.

### Usage

```python
from bhv_jpype import ensure_jvm, compute_bhv_distance, shutdown_jvm

ensure_jvm()
d = compute_bhv_distance("(A,(B,C));", "((A,B),C);", rooted=True)
print(d)
shutdown_jvm()
```
