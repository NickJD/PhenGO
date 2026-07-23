import os
import time
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


TARGET_FILES = [
    "MGI_GenePheno.rpt",
    "MGI_PhenoGenoMP.rpt",
    "VOC_MammalianPhenotype.rpt",
    "MGI_GenePhenoDisease.rpt",
]

# Older files were genuinely smaller — scale threshold by year
def min_expected_size_mb(year):
    if year <= 2012:
        return 0.1
    elif year <= 2015:
        return 0.5
    elif year <= 2018:
        return 1.0
    else:
        return 5.0


BAD_MARKERS = [
    "<html",
    "wayback machine",
    "internet archive",
    "blocked",
    "redirecting",
]


def make_session():
    """Session with automatic retries on transient network errors."""
    session = requests.Session()
    retry = Retry(
        total=5,
        backoff_factor=2,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "User-Agent": (
            "Mozilla/5.0 (X11; Linux x86_64) Chrome/120 Safari/537.36"
        )
    })
    return session


def is_valid_file(filepath, year):
    if not os.path.exists(filepath):
        return False

    size_mb = os.path.getsize(filepath) / (1024 * 1024)
    threshold = min_expected_size_mb(year)

    if size_mb < threshold:
        print(f"   Size {size_mb:.2f} MB below threshold {threshold} MB for {year}")
        return False

    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            head = f.read(2000).lower()
        for marker in BAD_MARKERS:
            if marker in head:
                print(f"   Rejected: found '{marker}' in content")
                return False
    except Exception:
        return False

    return True


def download_snapshot(session, raw_url, output_path):
    r = session.get(raw_url, stream=True, timeout=300, allow_redirects=True)
    r.raise_for_status()
    with open(output_path, "wb") as out:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                out.write(chunk)


def fetch_snapshots_cdx(session, target_file, status_filter=True):
    """
    Queries the Wayback CDX API across modern and historical legacy URL patterns.
    """
    cdx_url  = "https://web.archive.org/cdx/search/cdx"
    all_combined_snapshots = []

    # Map target file to all historical legacy alternative paths/names used by MGI
    url_variations = [
        f"http://www.informatics.jax.org/downloads/reports/{target_file}",
        f"https://jax.org{target_file}",
        f"ftp://ftp.informatics.jax.org/pub/reports/{target_file}"
    ]

    # Account for legacy name structural changes prior to 2015
    if target_file == "MGI_GenePhenoDisease.rpt":
        url_variations.append("ftp://ftp.informatics.jax.org/pub/reports/MGI_DiseasePheno.rpt")
        url_variations.append("http://jax.org")

    # De-duplicate variations to keep the network requests efficient
    url_variations = list(set(url_variations))

    for target_url in url_variations:
        params = {
            "url":       target_url,
            "output":    "json",
            "fl":        "timestamp,statuscode,original,length",
            "matchType": "exact",
            "limit":     "5000",
        }
        if status_filter:
            params["filter"] = "statuscode:200"

        try:
            print(f"  CDX query (status_filter={status_filter}): {target_url}")
            r = session.get(cdx_url, params=params, timeout=120)
            r.raise_for_status()
            data = r.json()
            
            if len(data) > 1:
                all_combined_snapshots.extend(data[1:])
        except Exception as e:
            print(f"  ⚠️ Variation check failed for {target_url}: {e}")
            continue
            
        time.sleep(1) # Polite gap between index checks

    return all_combined_snapshots


def organise_by_year(snapshots, start_year):
    """
    Group snapshots by year, keeping only those >= start_year.
    For each year, sort candidates: prefer mid-year (June) snapshots
    first, then fall back by descending size.
    """
    yearly = {}

    for row in snapshots:
        try:
            timestamp = row[0]
            original  = row[2]
            length    = int(row[3]) if len(row) > 3 else 0
            year      = int(timestamp[:4])
            month     = int(timestamp[4:6])

            if year < start_year:
                continue

            yearly.setdefault(year, []).append((timestamp, original, length, month))
        except Exception:
            continue

    # Sort: prefer snapshots closest to June (month 6), break ties by size desc
    for year in yearly:
        yearly[year].sort(key=lambda x: (abs(x[3] - 6), -x[2]))

    return yearly


def download_year(session, year, candidates, output_path):
    """Try each candidate snapshot until one succeeds."""
    for timestamp, original, length, month in candidates:
        raw_url     = f"https://web.archive.org/web/{timestamp}id_/{original}"
        expected_mb = round(length / (1024 * 1024), 1)

        print(f"\n  [{year}] Snapshot {timestamp} (month {month:02d}, "
              f"expected {expected_mb} MB)")
        print(f"  URL: {raw_url}")

        try:
            download_snapshot(session, raw_url, output_path)
            actual_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1)
            print(f"  Downloaded: {actual_mb} MB")

            if is_valid_file(output_path, year):
                print(f"  -> SUCCESS")
                return True
            else:
                print(f"  -> Invalid/truncated, trying next candidate")
                os.remove(output_path)

        except Exception as e:
            print(f"  -> Failed: {e}")
            if os.path.exists(output_path):
                os.remove(output_path)

        time.sleep(2)

    return False


def main(start_year=2010):
    session = make_session()

    for target_file in TARGET_FILES:
        file_base  = os.path.splitext(target_file)[0]
        output_dir = f"historical_mgi_{file_base}"
        os.makedirs(output_dir, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Processing: {target_file}")
        print(f"{'='*60}")

        # Primary CDX query (status 200 only)
        snapshots = fetch_snapshots_cdx(session, target_file, status_filter=True)
        print(f"Found {len(snapshots)} snapshots (status:200)")

        yearly = organise_by_year(snapshots, start_year)

        # Fallback: relaxed CDX query to fill any missing years
        if len(yearly) < (2025 - start_year):
            print("\nRunning relaxed CDX query to find missing years...")
            all_snapshots = fetch_snapshots_cdx(session, target_file, status_filter=False)
            all_yearly    = organise_by_year(all_snapshots, start_year)
            for year, candidates in all_yearly.items():
                if year not in yearly:
                    print(f"  Found fallback candidates for {year}")
                    yearly[year] = candidates

        covered   = sorted(yearly.keys())
        missing   = [y for y in range(start_year, 2026) if y not in yearly]
        print(f"\nYears with candidates : {covered}")
        if missing:
            print(f"Years with NO snapshots: {missing}")

        # Download
        results = {"success": [], "failed": []}

        for year in sorted(yearly):
            output_name = f"{file_base}_{year}.rpt"
            output_path = os.path.join(output_dir, output_name)

            # Skip already-valid downloads
            if is_valid_file(output_path, year):
                actual_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1)
                print(f"\n[{year}] Already have valid file ({actual_mb} MB) — skipping")
                results["success"].append(year)
                continue

            ok = download_year(session, year, yearly[year], output_path)
            (results["success"] if ok else results["failed"]).append(year)

            time.sleep(1)

        # Summary
        print(f"\n{'='*60}")
        print(f"Summary for {target_file}:")
        print(f"  Succeeded : {results['success']}")
        print(f"  Failed    : {results['failed']}")
        if missing:
            print(f"  No snapshots available: {missing}")
        print(f"{'='*60}\n")


if __name__ == "__main__":
    main(start_year=2010)

