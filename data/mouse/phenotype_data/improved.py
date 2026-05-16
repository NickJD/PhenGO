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

# Size thresholds scaled by year — older files were genuinely smaller
def min_expected_size_mb(year):
    if year <= 2012:
        return 0.05
    elif year <= 2015:
        return 0.2
    elif year <= 2018:
        return 0.5
    else:
        return 2.0


# HTML / stub markers that indicate a non-data response
BAD_CONTENT_MARKERS = [
    "<html",
    "<HTML",
    "wayback machine",
    "internet archive",
    "blocked",
    "redirecting",
    "not in archive",
    "no captures",
    "this page has not been",
    "404 not found",
    "access denied",
    "service unavailable",
    "{",          # JSON error body
    "<!doctype",
]

# MGI FTP — often has older files not captured by Wayback
FTP_BASE = "https://ftp.informatics.jax.org/pub/reports"

# All known URL variants for each file across MGI's history
URL_VARIANTS = [
    "http://www.informatics.jax.org/downloads/reports/{f}",
    "https://www.informatics.jax.org/downloads/reports/{f}",
    "http://informatics.jax.org/downloads/reports/{f}",
    "https://informatics.jax.org/downloads/reports/{f}",
]


def is_mgi_format(filepath):
    """Return True if the file looks like a real tab-separated MGI report."""
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for _ in range(20):
                line = f.readline()
                if not line:
                    break
                stripped = line.strip()
                if stripped.startswith("#") or stripped == "":
                    continue
                # A real MGI file will have multiple tab-separated columns
                if "\t" in stripped:
                    return True
                # Non-comment, non-empty line without tabs → not MGI format
                return False
    except Exception:
        pass
    return False


def is_valid_file(filepath, year):
    if not os.path.exists(filepath):
        return False

    size_mb = os.path.getsize(filepath) / (1024 * 1024)
    threshold = min_expected_size_mb(year)

    if size_mb < threshold:
        print(f"   Size {size_mb:.3f} MB below threshold {threshold} MB")
        return False

    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            head = f.read(4000)
        head_lower = head.lower()
        for marker in BAD_CONTENT_MARKERS:
            if marker.lower() in head_lower:
                print(f"   Rejected: stub marker '{marker}'")
                return False
    except Exception:
        return False

    if not is_mgi_format(filepath):
        print(f"   Rejected: not tab-separated MGI data")
        return False

    return True


def make_session():
    session = requests.Session()
    retry = Retry(
        total=6,
        backoff_factor=3,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Chrome/120 Safari/537.36"
    })
    return session


def download_snapshot(session, raw_url, output_path):
    r = session.get(raw_url, stream=True, timeout=300, allow_redirects=True)
    r.raise_for_status()
    with open(output_path, "wb") as out:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                out.write(chunk)


def fetch_cdx_for_url(session, url, status_filter):
    params = {
        "url":       url,
        "output":    "json",
        "fl":        "timestamp,statuscode,original,length",
        "matchType": "exact",
        "limit":     "5000",
    }
    if status_filter:
        params["filter"] = "statuscode:200"

    r = session.get(
        "https://web.archive.org/cdx/search/cdx",
        params=params,
        timeout=120,
    )
    r.raise_for_status()
    data = r.json()
    return data[1:] if len(data) > 1 else []


def fetch_all_snapshots(session, target_file):
    """
    Query CDX across all known URL variants, with and without status filter.
    De-duplicate by timestamp.
    """
    seen     = set()
    all_rows = []

    for variant_tpl in URL_VARIANTS:
        url = variant_tpl.format(f=target_file)
        for status_filter in (True, False):
            try:
                rows = fetch_cdx_for_url(session, url, status_filter)
                for row in rows:
                    ts = row[0]
                    if ts not in seen:
                        seen.add(ts)
                        all_rows.append(row)
            except Exception as e:
                print(f"   CDX failed for {url} (filter={status_filter}): {e}")

    print(f"  Unique snapshots across all URL variants: {len(all_rows)}")
    return all_rows


def organise_by_year(snapshots, start_year):
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

    # Prefer mid-year (June) snapshots; break ties by descending size
    for year in yearly:
        yearly[year].sort(key=lambda x: (abs(x[3] - 6), -x[2]))

    return yearly


def try_ftp_direct(session, target_file, output_path, year):
    """Fetch directly from MGI's FTP HTTPS gateway (useful for current/recent files)."""
    url = f"{FTP_BASE}/{target_file}"
    print(f"\n  Trying FTP direct: {url}")
    try:
        download_snapshot(session, url, output_path)
        actual_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1)
        print(f"  Downloaded: {actual_mb} MB")
        if is_valid_file(output_path, year):
            print(f"  -> FTP SUCCESS")
            return True
        os.remove(output_path)
    except Exception as e:
        print(f"  -> FTP failed: {e}")
        if os.path.exists(output_path):
            os.remove(output_path)
    return False


def download_year(session, year, candidates, output_path):
    """Try each Wayback candidate snapshot until one succeeds."""
    for timestamp, original, length, month in candidates:
        raw_url     = f"https://web.archive.org/web/{timestamp}id_/{original}"
        expected_mb = round(length / (1024 * 1024), 1)

        print(f"\n  [{year}] Snapshot {timestamp} "
              f"(month {month:02d}, expected {expected_mb} MB)")

        try:
            download_snapshot(session, raw_url, output_path)
            actual_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1)
            print(f"  Downloaded: {actual_mb} MB")

            if is_valid_file(output_path, year):
                print(f"  -> SUCCESS")
                return True
            else:
                print(f"  -> Rejected, trying next")
                os.remove(output_path)

        except Exception as e:
            print(f"  -> Failed: {e}")
            if os.path.exists(output_path):
                os.remove(output_path)

        time.sleep(2)

    return False


def main(start_year=2010):
    session      = make_session()
    current_year = 2026

    for target_file in TARGET_FILES:
        file_base  = os.path.splitext(target_file)[0]
        output_dir = f"historical_mgi_{file_base}"
        os.makedirs(output_dir, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Processing: {target_file}")
        print(f"{'='*60}")

        snapshots = fetch_all_snapshots(session, target_file)
        yearly    = organise_by_year(snapshots, start_year)

        covered = sorted(yearly.keys())
        missing = [y for y in range(start_year, current_year + 1) if y not in yearly]
        print(f"Years with candidates  : {covered}")
        print(f"Years with no snapshots: {missing}")

        results = {
            "success":     [],
            "failed":      [],
            "no_snapshot": list(missing),
        }

        for year in sorted(yearly):
            output_name = f"{file_base}_{year}.rpt"
            output_path = os.path.join(output_dir, output_name)

            if is_valid_file(output_path, year):
                actual_mb = round(os.path.getsize(output_path) / (1024 * 1024), 1)
                print(f"\n[{year}] Already valid ({actual_mb} MB) — skipping")
                results["success"].append(year)
                continue

            ok = download_year(session, year, yearly[year], output_path)
            (results["success"] if ok else results["failed"]).append(year)
            time.sleep(1)

        # FTP fallback for current year if Wayback failed or had no snapshot
        if current_year in results["failed"] or current_year in results["no_snapshot"]:
            output_path = os.path.join(
                output_dir, f"{file_base}_{current_year}.rpt"
            )
            if not is_valid_file(output_path, current_year):
                ok = try_ftp_direct(session, target_file, output_path, current_year)
                if ok:
                    for bucket in ("failed", "no_snapshot"):
                        if current_year in results[bucket]:
                            results[bucket].remove(current_year)
                    results["success"].append(current_year)

        print(f"\n{'='*60}")
        print(f"Summary for {target_file}:")
        print(f"  Succeeded        : {sorted(results['success'])}")
        print(f"  Failed (tried)   : {sorted(results['failed'])}")
        print(f"  No archive found : {sorted(results['no_snapshot'])}")
        print(f"{'='*60}")

    print(
        "\nNote: 'No archive found' years were never captured by the Wayback Machine.\n"
        "For pre-2016 MGI data, try contacting MGI directly or browsing:\n"
        f"  {FTP_BASE}/"
    )


if __name__ == "__main__":
    main(start_year=2010)
