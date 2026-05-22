#!/usr/bin/env python3

import requests
import json
import os

def D1_MGnify_search(max_pages, url=None, output_prefix="aquatic", output_dir="MGnify_pipeline_results"):
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

##############################################################################################

def D2_genome_urls_list(id_list, output_file="metagenomes_url_links.txt", output_dir = "MGnify_pipeline_results"):
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

##############################################################################################

def D3_download_genome_url(input, filter_ext=".fna"):
    """
    Reads a .txt file OR a list of links and downloads only the desired files by filtering by file type.
    
    Parameters:
    input: path to a .txt file containing the links OR a list of links.
    filter_ext: file extension to be downloaded. Default: ".fna"
    """

    output_dir = "./metagenomes_MGnify"
    os.makedirs(output_dir, exist_ok=True)

    links = [line.strip() for line in input if line.strip()]
    filtered_links = [link for link in links if link.endswith(filter_ext)] 

    if not filtered_links:
        print(f"No {filter_ext} links found.") 
        return

    for link in filtered_links:
        filename = link.split("/")[-1]
        print(filename)
        print(f"Downloading {filename} ...")

        try:
            response = requests.get(link)
            response.raise_for_status()
            with open(f"./metagenomes_MGnify/{filename}", "wb") as f:
                f.write(response.content)
            print(f"Download completed: {filename}")
        except requests.RequestException as e:
            print(f"Error downloading {filename}: {e}")