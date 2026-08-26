"""Driver audit for the conversion-risk SI analysis.

The conversion hazard used in `si_deforestation.R` is GFW tree-cover loss restricted
to five non-fire anthropogenic drivers under the WRI/Google DeepMind 1 km dominant-driver
attribution of Sims et al. (2025).
Loss attributed to Wildfire, Other natural disturbances or Unknown is excluded. This
script measures how much loss falls into each bucket, for two reasons:

  1. `Unknown` is loss that enters neither the fire hazard nor the conversion hazard,
     so it biases both downwards. The SI must state its size.
  2. The `Wildfire` share quantifies the double-counting that the driver attribution
     avoids. Using undifferentiated tree-cover loss alongside the EFFIS fire model
     would count that share twice.

Queries the GFW Data API's public /download/json endpoint. No authentication.

Self-contained: it imports nothing from outside this repository and writes only to
`outputs/si/`.

Query parameters are pinned to match the panel the simulation reads
(`data/conversion_rates.csv`), so the shares reported here describe that exact file.

Run from repo root:
    python3 code/si/si_defor_driver_audit.py
"""
from __future__ import annotations

import json
import sys
import urllib.parse
from pathlib import Path

import pandas as pd
import requests

# These are the values the shipped conversion panel was built with; the audit is only
# meaningful if it queries the same dataset version, threshold and year range.
BASE_URL = "https://data-api.globalforestwatch.org"
CHANGE_DATASET = "gadm__tcl__adm1_change"
VERSION = "v20260424"
CANOPY_THRESHOLD = 30
START_YEAR, END_YEAR = 2001, 2025

NON_FIRE_DRIVERS = [
    "Hard commodities",
    "Logging",
    "Permanent agriculture",
    "Settlements & Infrastructure",
    "Shifting cultivation",
]
FIRE_DRIVER = "Wildfire"
NATURAL_DRIVER = "Other natural disturbances"
UNKNOWN_DRIVER = "Unknown"

# (iso, adm1) as used in jurisdiction_annual_rates.csv, matching the EFFIS cache keys
REGIONS = {("USA", 5): "California", ("BRA", 12): "Mato Grosso", ("IDN", 23): "Papua"}

OUT_DIR = Path("outputs/si")
CACHE_DIR = Path("data/gfw_cache")


def fetch_change(iso: str) -> pd.DataFrame:
    """Annual tree-cover loss by adm1, year and driver for one country."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / f"{iso}_change.json"

    if cache.exists():
        rows = json.loads(cache.read_text())
        print(f"  [{iso}] {len(rows)} rows from cache")
    else:
        sql = (
            "SELECT iso, adm1, umd_tree_cover_loss__year AS year, "
            "wri_google_tree_cover_loss_drivers__driver AS driver, "
            "SUM(umd_tree_cover_loss__ha) AS loss_ha "
            f"FROM data WHERE iso = '{iso}' "
            f"AND umd_tree_cover_density_2000__threshold = {CANOPY_THRESHOLD} "
            "GROUP BY iso, adm1, umd_tree_cover_loss__year, "
            "wri_google_tree_cover_loss_drivers__driver"
        )
        url = (f"{BASE_URL}/dataset/{CHANGE_DATASET}/{VERSION}/download/json"
               f"?sql={urllib.parse.quote(sql)}")
        print(f"  [{iso}] fetching ...")
        resp = requests.get(url, timeout=180)
        resp.raise_for_status()
        rows = resp.json()
        cache.write_text(json.dumps(rows))
        print(f"  [{iso}] {len(rows)} rows -> {cache}")

    df = pd.DataFrame(rows)
    df["loss_ha"] = pd.to_numeric(df["loss_ha"], errors="coerce")
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    df["adm1"] = pd.to_numeric(df["adm1"], errors="coerce")
    return df[(df.year >= START_YEAR) & (df.year <= END_YEAR)]


def bucket(driver: str) -> str:
    if driver in NON_FIRE_DRIVERS:
        return "non_fire_anthropogenic"
    if driver == FIRE_DRIVER:
        return "wildfire"
    if driver == NATURAL_DRIVER:
        return "other_natural"
    return "unknown"


def main() -> int:
    print("GFW driver audit for the three focal regions")
    print(f"dataset {CHANGE_DATASET} {VERSION}, canopy {CANOPY_THRESHOLD}%, "
          f"{START_YEAR}-{END_YEAR}\n")

    frames = {iso: fetch_change(iso) for iso in sorted({i for i, _ in REGIONS})}

    records = []
    for (iso, adm1), name in REGIONS.items():
        sub = frames[iso]
        sub = sub[sub.adm1 == adm1].copy()
        if sub.empty:
            print(f"\n!! no rows for {name} ({iso} adm1={adm1})")
            continue
        sub["bucket"] = sub.driver.map(bucket)
        tot = sub.loss_ha.sum()
        by = sub.groupby("bucket").loss_ha.sum()
        rec = {"region": name, "iso": iso, "adm1": adm1, "total_loss_ha": tot}
        for b in ["non_fire_anthropogenic", "wildfire", "other_natural", "unknown"]:
            rec[f"{b}_ha"] = float(by.get(b, 0.0))
            rec[f"{b}_share"] = float(by.get(b, 0.0)) / tot if tot else float("nan")
        records.append(rec)

        # Drivers actually present, to catch a renamed category upstream
        unmapped = sorted(set(sub[sub.bucket == "unknown"].driver.unique()) - {UNKNOWN_DRIVER})
        if unmapped:
            print(f"\n!! {name}: drivers falling into 'unknown' that are not "
                  f"'{UNKNOWN_DRIVER}': {unmapped}")

    out = pd.DataFrame(records)

    print("\nShare of total tree-cover loss, 2001-2025:\n")
    hdr = f"{'region':<13}{'non-fire (delta)':>18}{'wildfire':>11}{'natural':>10}{'unknown':>10}"
    print(hdr)
    print("-" * len(hdr))
    for _, r in out.iterrows():
        print(f"{r.region:<13}{100*r.non_fire_anthropogenic_share:>17.2f}%"
              f"{100*r.wildfire_share:>10.2f}%{100*r.other_natural_share:>9.2f}%"
              f"{100*r.unknown_share:>9.2f}%")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    dest = OUT_DIR / "si_defor_driver_audit.csv"
    out.to_csv(dest, index=False)
    print(f"\nWrote {dest}")

    print("\nFor the SI:")
    for _, r in out.iterrows():
        print(f"  {r.region}: wildfire is {100*r.wildfire_share:.1f}% of tree-cover loss "
              f"(would be double-counted without driver attribution); "
              f"unknown is {100*r.unknown_share:.2f}% (excluded from both hazards).")
    worst = out.unknown_share.max()
    print(f"\nLargest unknown share: {100*worst:.2f}% -> "
          f"{'state as a quantified downward bias' if worst > 0.02 else 'immaterial, note in passing'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
