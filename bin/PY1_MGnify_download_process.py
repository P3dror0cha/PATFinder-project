#!/usr/bin/env python3

import argparse

from PY1_MGnify_download_functions import D1_MGnify_search
from PY1_MGnify_download_functions import D2_genome_urls_list
from PY1_MGnify_download_functions import D3_download_genome_url

def MGnify_pipeline(max_pages, output_prefix="aquatic", extension=".fna", output_dir="./genomes_MGnify"):
    """
    The complete MGnify pipeline integrates the following functions:
    1. D1_MGnify_search
    2. D2_genome_urls_list
    3. D3_download_genome_url

    Parameters:
    max_pages (int): Number of pages to retrieve from MGnify. Each page has 25 metagenomes. (required)
    output_prefix (str): Prefix for intermediate output files.
    extension (str): File extension to filter downloads.
    output_dir (str): Directory to save downloaded files.

    Returns:
    None
    """

    json_file, id_list = D1_MGnify_search(max_pages)
    url_list = D2_genome_urls_list(id_list, output_file=f"{output_prefix}_url_list.txt")
    D3_download_genome_url(url_list, filter_ext=".fna")
    D3_download_genome_url(url_list, filter_ext=".faa")



def main():
    parser = argparse.ArgumentParser(
        description="Complete MGnify genome download pipeline"
    )

    parser.add_argument(
        "--max_pages",
        type=int,
        required=True,
        help="Number of pages to retrieve from MGnify (required)"
    )

    parser.add_argument(
        "--output_prefix",
        default="aquatic",
        help="Prefix for intermediate output files"
    )

    args = parser.parse_args()

    MGnify_pipeline(
        max_pages=args.max_pages,
        output_prefix=args.output_prefix
    )

if __name__ == "__main__":
    main()