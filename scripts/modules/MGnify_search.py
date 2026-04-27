#!/usr/bin/env python3

import requests
import json
import argparse
import os

def MGnify_search(max_pages, url=None, output_prefix="aquatic", output_dir="MGnify_pipeline_results"):
    """
    Retrieves multiple sample IDs from the MGnify database.

    Parameters:
    max_pages: maximum number of pages to retrieve.
    url: custom MGnify API URL (optional).
    output_prefix: prefix for output files.
    output_dir: Output directory for the id list and the json file.
    """

    if url is None:
        url = "https://www.ebi.ac.uk/metagenomics/api/v1/biomes/root:Environmental:Aquatic:Marine/genomes"

    print("Starting ID download")

    json_list = []
    id_list = []

    os.makedirs(output_dir, exist_ok=True)

    for page in range(1, max_pages + 1):
        response = requests.get(url, params={"page": page})

        if response.status_code != 200:
            print(f"Error {response.status_code} on page {page}")
            break

        data = response.json().get("data", [])

        json_list.extend(data)
        id_list.extend(item["id"] for item in data)

        print(f"Page {page}: {len(data)} IDs found")

    final_json = {"data": json_list}

    json_file = f"{output_prefix}_download.json"
    ids_file = f"{output_prefix}_ids.txt"

    with open(f"./MGnify_pipeline_results/{json_file}", "w") as f_json, open(f"./MGnify_pipeline_results/{ids_file}", "w") as f_ids:
        json.dump(final_json, f_json, indent=4)
        f_ids.write("\n".join(id_list))

    print(f"\nTotal of {len(id_list)} IDs collected.")
    print("End of ID download")

    return final_json, id_list


def main():
    parser = argparse.ArgumentParser(
        description="Download IDs from MGnify"
    )

    parser.add_argument(
        "--max_pages",
        type=int,
        required=True,
        help="Maximum number of pages to retrieve"
    )

    parser.add_argument(
        "--url",
        default=None,
        help="Custom MGnify API URL"
    )

    parser.add_argument(
        "--output_prefix",
        default="aquatic",
        help="Prefix for output files"
    )

    parser.add_argument(
    "--output_dir",
    default="MGnify_pipeline_results",
    help="Name for the output directory."
    )

    args = parser.parse_args()

    MGnify_search(
        max_pages=args.max_pages,
        url=args.url,
        output_prefix=args.output_prefix,
        output_dir=args.output_dir 
    )


if __name__ == "__main__":
    main()
