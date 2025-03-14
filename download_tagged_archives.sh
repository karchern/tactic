#!/bin/bash

# Define the repository and tags
REPO="karchern/tactic"
TAGS=("screen_1_WITHOUT_wrong_clones" "screen_1_WITH_wrong_clones" "screen_2")  # Add your tags here

# Loop through each tag and download the archive
for TAG in "${TAGS[@]}"; do
  URL="https://github.com/$REPO/archive/refs/tags/$TAG.tar.gz"
  wget -O "tagged_archives/$TAG.tar.gz" "$URL"
  # unpack
  mkdir -p tagged_archives_unpacked
  tar -xvf "tagged_archives/$TAG.tar.gz" -C "tagged_archives_unpacked/"
done