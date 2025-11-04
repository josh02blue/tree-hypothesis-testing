import pandas as pd
import os

# === Configuration ===
CSV_PATH = "./batch_output/new_permutation_summary_304.csv"   # path to your combined CSV
ALPHA = 0.01      #1e-3                         # p-value threshold
EPS = 0.10        #0.25                            # margin as % of (null_max - null_min)
ABS_TOL = 0.0                          # minimum absolute margin
OUTPUT_FILE = "new_odd_cases.csv"          # where to save results

# === Load CSV ===
df = pd.read_csv(CSV_PATH)

# --- Map columns to lowercase names if needed ---
rename_map = {
    "Locus": "filename",
    "Test_Statistic": "test_stat",
    "Null_Min": "null_min",
    "Null_Max": "null_max",
    "Null_Mean": "null_mean",
    "P_Value": "p_value",
}
df = df.rename(columns={c: rename_map[c] for c in rename_map if c in df.columns})

# Normalize filename column
if "filename" in df.columns:
    df["filename"] = df["filename"].astype(str).apply(lambda x: os.path.basename(str(x)).strip())

# === Column names ===
name = "filename"
test = "test_stat"
pcol = "p_value"
nmax = "null_max"
nmin = "null_min"

# --- Outside-null rule with margin ---
if nmax in df.columns and nmin in df.columns:
    rng = (df[nmax] - df[nmin]).abs()
    margin = (EPS * rng).clip(lower=ABS_TOL)
    outside_hi = df[test] >= (df[nmax] + margin)
    outside_lo = df[test] <= (df[nmin] - margin)
    outside = outside_hi | outside_lo
else:
    margin = pd.Series([float("nan")] * len(df))
    outside = pd.Series([False] * len(df))

# --- Significance rule ---
sig = df[pcol] < ALPHA if pcol in df.columns else pd.Series([True] * len(df))

# --- Final mask for odd cases ---
mask = outside & sig
hits = df.loc[mask, name].astype(str).tolist()

# --- Summary table ---
summary_cols = [c for c in [name, test, nmin, nmax, pcol] if c in df.columns]
summary = df.loc[mask, summary_cols].copy()

if nmax in df.columns and nmin in df.columns:
    summary["margin_used"] = margin[mask].values
    summary["outside_null"] = True
summary["p_sig"] = True

# --- Output ---
print("\n✅ Odd cases (filenames):")
for h in hits:
    print(h)

summary.to_csv(OUTPUT_FILE, index=False)
print(f"\n📁 Saved table → {OUTPUT_FILE}")
print(f"Total odd cases found: {len(summary)}")