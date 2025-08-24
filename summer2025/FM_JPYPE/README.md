# Fréchet Mean via JPype (SturmMean_201102.jar)

A lightweight Python wrapper that calls the Fréchet mean algorithm in
**SturmMean_201102.jar** using **JPype**. Works standalone or as part of
a pipeline with BHV (gtp.jar).

## Files

- `frechet_jpype.py` — importable module (start/stop JVM + `frechet_mean_jpype`)
- `example_frechet.py` — tiny example that prints a mean tree

## Requirements

- Python ≥ 3.8
- Java (JRE/JDK) available
- `jpype1` (`pip install jpype1`)
- `SturmMean_201102.jar` available in this folder (or update the path in the module)

## Quick start

```bash
pip install jpype1
python example_frechet.py
```
