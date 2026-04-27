#!/usr/bin/env python3

import requests
import argparse

def genome_urls_list(id_list, output_file="metagenomes_url_links.txt", output_dir = "MGnify_pipeline_results"):
    """
    Retrieve download URLs associated with each genome ID from the MGnify API
    and save them into a single text file.
    """

    download_links = []

    print("Starting download URL retrieval...\n")

    for genome_id in id_list:
        print(f"Fetching links for ID: {genome_id}")

        url = f"https://www.ebi.ac.uk/metagenomics/api/v1/genomes/{genome_id}/downloads"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            data = response.json().get("data", [])
            urls = [item["links"]["self"] for item in data]

            download_links.extend(urls)

            print(f"{len(urls)} links found\n")

        except requests.exceptions.RequestException as e:
            print(f"Request failed for ID {genome_id}: {e}")

    with open(f"./{output_dir}/{output_file}", "w") as f:
        f.write("\n".join(download_links))

    print(f"\nTotal links collected: {len(download_links)}")
    print(f"Saved to: {output_file}")

    return download_links


def main():
    parser = argparse.ArgumentParser(
        description="Fetch download URLs from MGnify genome IDs"
    )

    parser.add_argument(
        "--input_ids",
        required=True,
        help="File containing genome IDs (one per line)"
    )

    parser.add_argument(
        "--output_file",
        default="metagenome_download_links.txt",
        help="Output file name"
    )

    parser.add_argument(
        "--output_dir",
        default="MGnify_pipeline_results",
        help="Output file name"
    )
    args = parser.parse_args()

    with open(args.input_ids) as f:
        id_list = [line.strip() for line in f if line.strip()]

    genome_urls_list(
        id_list=id_list,
        output_file=args.output_file,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()