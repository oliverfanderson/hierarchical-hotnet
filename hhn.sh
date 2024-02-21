#!/usr/bin/env bash

# Initialize variables
network=""
scores=""
data=$PWD/data
intermediate=$PWD/intermediate
results=$PWD/results

num_permutations=100

# Parse command-line options
while getopts ":n:s:" opt; do
  case $opt in
    n) network="$OPTARG" ;;
    s) scores="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# Check if required options are provided
if [ -z "$network" ] || [ -z "$scores" ]; then
  echo "Usage: $0 -n <network> -s <scores>"
  exit 1
fi

# Convert comma-separated scores to array
IFS=',' read -r -a scores_array <<< "$scores"

# Hierarchical HotNet is parallelizable, but this script runs each Hierarchical
# HotNet sequentially.  Please see the example_commands_parallel.sh script for a
# parallelized example.

# Compile Fortran module.
cd ./src
f2py -c fortran_module.f95 -m fortran_module > /dev/null
cd ..

################################################################################
#
#   Prepare data.
#
################################################################################

# Create data, intermediate data and results, and results directories.
mkdir --p "$data"
mkdir --p "$intermediate"
mkdir --p "$results"

# Create directories for network
mkdir --p "$intermediate/$network"

# Create directories for network and scores
for score in "${scores_array[@]}"; do
  mkdir --p "$intermediate/"$network"_"$score""
done

################################################################################
#
#   Construct similarity matrices.
#
################################################################################

for network in "$network"; do
  python src/construct_similarity_matrix.py \
    -i   "$data/"$network"_edge_list.tsv" \
    -o   "$intermediate/"$network"/similarity_matrix.h5" \
    -bof "$intermediate/"$network"/beta.txt"
done

################################################################################
#
#   Permute data.
#
################################################################################

# This example does not use permuted networks, but these commands show how to
# generate them.
echo "Permuting networks..."

for network in "$network"; do
  for score in "${scores_array[@]}"; do
    cp "$data/"$network"_index_gene.tsv" "$intermediate/"$network"/index_gene_0.tsv"
    cp "$data/"$network"_edge_list.tsv" "$intermediate/"$network"/edge_list_0.tsv"

    # Preserve connectivity of the observed graph.
    for i in `seq 1 4`; do
      python src/permute_network.py \
        -i "$intermediate/"$network"/edge_list_0.tsv" \
        -s "$i" \
        -c \
        -o "$intermediate/"$network"/edge_list_${i}.tsv"
    done

    # Do not preserve connectivity of the observed graph.
    for i in `seq 5 8`; do
      python src/permute_network.py \
        -i "$intermediate/"$network"/edge_list_0.tsv" \
        -s "$i" \
        -o "$intermediate/"$network"/edge_list_${i}.tsv"
    done
  done
done

echo "Permuting scores..."

for network in "$network"; do
  for score in "${scores_array[@]}"; do
    cp "$data/"$score".tsv" "$intermediate/"$network"_"$score"/scores_0.tsv"

    python src/find_permutation_bins.py \
      -gsf "$intermediate/"$network"_"$score"/scores_0.tsv" \
      -igf "$data/"$network"_index_gene.tsv" \
      -elf "$data/"$network"_edge_list.tsv" \
      -ms 1000 \
      -o "$intermediate/"$network"_"$score"/score_bins.tsv"

    for i in `seq $num_permutations`; do
      python src/permute_scores.py \
        -i "$intermediate/"$network"_"$score"/scores_0.tsv" \
        -bf "$intermediate/"$network"_"$score"/score_bins.tsv" \
        -s "$i" \
        -o "$intermediate/"$network"_"$score"/scores_${i}.tsv"
    done
  done
done

################################################################################
#
#   Construct hierarchies.
#
################################################################################

echo "Constructing hierarchies..."

for network in "$network"; do
  for score in "${scores_array[@]}"; do
    for i in $(seq 0 $num_permutations); do
      python src/construct_hierarchy.py \
        -smf "$intermediate/"$network"/similarity_matrix.h5" \
        -igf "$data/"$network"_index_gene.tsv" \
        -gsf "$intermediate/"$network"_"$score"/scores_${i}.tsv" \
        -helf "$intermediate/"$network"_"$score"/hierarchy_edge_list_${i}.tsv" \
        -higf "$intermediate/"$network"_"$score"/hierarchy_index_gene_${i}.tsv"
    done
  done
done

################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

# This example uses -lsb/--lower_size_bound 1 because it is a small toy example
# with 25 vertices.  Use larger value (default is 10) for larger graphs.
for network in "$network"; do
  for score in "${scores_array[@]}"; do
    python src/process_hierarchies.py \
      -oelf "$intermediate/"$network"_"$score"/hierarchy_edge_list_0.tsv" \
      -oigf "$intermediate/"$network"_"$score"/hierarchy_index_gene_0.tsv" \
      -pelf $(for i in $(seq $num_permutations); do echo "$intermediate/"$network"_"$score"/hierarchy_edge_list_${i}.tsv "; done) \
      -pigf $(for i in $(seq $num_permutations); do echo "$intermediate/"$network"_"$score"/hierarchy_index_gene_${i}.tsv "; done) \
      -lsb 1 \
      -cf "$results/clusters_"$network"_"$score".tsv" \
      -pl "$network" "$score" \
      -pf "$results/sizes_"$network"_"$score".pdf"
  done
done

################################################################################
#
#   Perform consensus.
#
################################################################################

# echo "Performing consensus..."

# python src/perform_consensus.py \
#   -cf "$results/clusters_"$network"_"${scores_array[0]}".tsv" "$results/clusters_"$network"_"${scores_array[1]}".tsv" \
#   -igf "$data/"$network"_index_gene.tsv" "$data/"$network"_index_gene.tsv" \
#   -elf "$data/"$network"_edge_list.tsv" "$data/"$network"_edge_list.tsv" \
#   -n "$network" "$network" \
#   -s "${scores_array[@]}" \
#   -t 2 \
#   -cnf "$results/consensus_nodes.tsv" \
#   -cef "$results/consensus_edges.tsv"