import os
import time
import requests


TARGET_FILES = [
    #"MGI_GenePheno.rpt",
    #"MGI_PhenoGenoMP.rpt",
    "VOC_MammalianPhenotype.rpt",
    "MGI_GenePhenoDisease.rpt"
]


MIN_EXPECTED_SIZE_MB = 5


def is_valid_file(filepath):

    if not os.path.exists(filepath):
        return False

    size_mb = os.path.getsize(filepath) / (1024 * 1024)

    if size_mb < MIN_EXPECTED_SIZE_MB:
        return False

    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            head = f.read(2000).lower()

        bad_markers = [
            "<html",
            "wayback machine",
            "internet archive",
            "blocked",
            "redirecting"
        ]

        for marker in bad_markers:
            if marker in head:
                return False

    except Exception:
        return False

    return True


def download_snapshot(raw_url, output_path, headers):

    r = requests.get(
        raw_url,
        headers=headers,
        stream=True,
        timeout=300,
        allow_redirects=True
    )

    r.raise_for_status()

    with open(output_path, "wb") as out:
        for chunk in r.iter_content(chunk_size=1024 * 1024):
            if chunk:
                out.write(chunk)


def fetch_snapshots_for_file(target_file):

    base_url = (
        f"http://www.informatics.jax.org/downloads/reports/{target_file}"
    )

    cdx_url = "https://web.archive.org/cdx/search/cdx"

    params = {
        "url": base_url,
        "output": "json",
        "fl": "timestamp,statuscode,original,length",
        "filter": "statuscode:200",
        "matchType": "exact"
    }

    headers = {
        "User-Agent": (
            "Mozilla/5.0 "
            "(X11; Linux x86_64) "
            "Chrome/120 Safari/537.36"
        )
    }

    print(f"\nSearching snapshots:")
    print(base_url)

    r = requests.get(
        cdx_url,
        params=params,
        headers=headers,
        timeout=120
    )

    r.raise_for_status()

    data = r.json()

    if len(data) <= 1:
        return []

    return data[1:]


def main(start_year=2010):

    headers = {
        "User-Agent": (
            "Mozilla/5.0 "
            "(X11; Linux x86_64) "
            "Chrome/120 Safari/537.36"
        )
    }

    for target_file in TARGET_FILES:

        file_base = os.path.splitext(target_file)[0]

        output_dir = f"historical_mgi_{file_base}"

        os.makedirs(output_dir, exist_ok=True)

        snapshots = fetch_snapshots_for_file(target_file)

        print(f"Found {len(snapshots)} snapshots")

        # Organise by year
        yearly_snapshots = {}

        for row in snapshots:

            try:
                timestamp = row[0]
                status = row[1]
                original = row[2]

                length = 0

                if len(row) > 3:
                    try:
                        length = int(row[3])
                    except:
                        pass

                year = int(timestamp[:4])

                if year < start_year:
                    continue

                yearly_snapshots.setdefault(year, []).append(
                    (timestamp, original, length)
                )

            except Exception:
                continue

        for year in sorted(yearly_snapshots):

            candidates = sorted(
                yearly_snapshots[year],
                key=lambda x: x[2],
                reverse=True
            )

            success = False

            for timestamp, original, length in candidates:

                output_name = f"{file_base}_{year}.rpt"

                output_path = os.path.join(
                    output_dir,
                    output_name
                )

                raw_url = (
                    f"https://web.archive.org/web/"
                    f"{timestamp}id_/"
                    f"{original}"
                )

                expected_mb = round(length / (1024 * 1024), 1)

                print(f"\n[{year}] Trying snapshot:")
                print(timestamp)
                print(f"Expected archive size: {expected_mb} MB")

                try:

                    download_snapshot(
                        raw_url,
                        output_path,
                        headers
                    )

                    actual_mb = round(
                        os.path.getsize(output_path)
                        / (1024 * 1024),
                        1
                    )

                    print(f"Downloaded size: {actual_mb} MB")

                    if not is_valid_file(output_path):

                        print("-> Invalid/truncated snapshot")

                        os.remove(output_path)

                        continue

                    print(f"-> SUCCESS")
                    success = True
                    break

                except Exception as e:

                    print(f"-> Failed: {e}")

                    if os.path.exists(output_path):
                        os.remove(output_path)

                time.sleep(2)

            if not success:
                print(f"[{year}] No valid snapshot found")


if __name__ == "__main__":
    main(start_year=2010)
