"""Download GSE161529 raw data from GEO."""

import logging
import sys

from trajectory.io import download_geo

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def main() -> None:
    force = "--force" in sys.argv
    dest = download_geo(force=force)
    print(f"Data available at: {dest}")


if __name__ == "__main__":
    main()
