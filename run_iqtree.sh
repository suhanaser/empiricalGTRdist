threads=26

# run IQ tree on every dataset to find the best-fit GTR matrix and base frequencies for each partition.
folders=$(find /data/Suha/GTR_parameters_dist/ -type f -name 'alignment.nex' -printf '%h\n' | sort -u)
echo "$folders" | parallel -P $threads iqtree -s {}"/alignment.nex" -spp {}"/alignment.nex"  -n 0 -m GTR+I --prefix {}"/model"

# run IQ tree on every dataset to find the ML topology and branch lengths.
folders=$(find /data/Suha/GTR_parameters_dist/ -type f -name 'alignment.nex' -printf '%h\n' | sort -u)
echo "$folders" | parallel -P $threads iqtree -s {}"/alignment.nex" -spp {}"/alignment.nex" --prefix {}"/branches"