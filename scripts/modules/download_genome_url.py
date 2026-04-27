#!/usr/bin/env python3

import requests
import os

def download_genome_url(input, filter_ext=".fna"):
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
