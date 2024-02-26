#!/usr/bin/env bash

# Initialize variables
scores=""
edges=""
index=""
output_file=""
# data=$PWD/data
intermediate="HHN/intermediate"
results="HHN/results"

num_permutations=100

# Parse command-line options
while getopts ":s:e:i:o:" opt; do
  case $opt in
    s) scores="$OPTARG" ;;
    e) edges="$OPTARG" ;;
    i) index="$OPTARG" ;;
    o) output_file="$OPTARG" ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# Check if required options are provided
if [ -z "$scores" ] || [ -z "$edges" ] || [ -z "$index" ] || [ -z "$output_file" ]; then
  echo "Usage: $0 -s <scores> -e <edges> -i <index> -o <output_file>"
  exit 1
fi

# Hierarchical HotNet is parallelizable, but this script runs each Hierarchical
# HotNet sequentially.  Please see the example_commands_parallel.sh script for a
# parallelized example.

# Compile Fortran module.
cd HHN/src
f2py -c fortran_module.f95 -m fortran_module > /dev/null
cd ..

################################################################################
#
#   Prepare data.
#
################################################################################

# Create data, intermediate data and results, and results directories.
# mkdir --p "$data"
mkdir -p "$intermediate"
mkdir -p "$results"

# Create directories for network
mkdir -p "$intermediate/network"

# Create directories for network and scores
mkdir -p "$intermediate/network_score"

################################################################################
#
#   Construct similarity matrices.
#
################################################################################


python $PWD/src/construct_similarity_matrix.py \
	-i   "$edges" \
	-o   "$intermediate/network/similarity_matrix.h5" \
	-bof "$intermediate/network/beta.txt"

################################################################################
#
#   Permute data.
#
################################################################################

# This example does not use permuted networks, but these commands show how to
# generate them.
echo "Permuting networks..."

cp "$index" "$intermediate/network/index_gene_0.tsv"
cp "$edges" "$intermediate/network/edge_list_0.tsv"

    # Preserve connectivity of the observed graph.
    for i in `seq 1 4`; do
      python $PWD/src/permute_network.py \
        -i "$intermediate/network/edge_list_0.tsv" \
        -s "$i" \
        -c \
        -o "$intermediate/network/edge_list_${i}.tsv"
    done

    # Do not preserve connectivity of the observed graph.
    for i in `seq 5 8`; do
      python $PWD/src/permute_network.py \
        -i "$intermediate/network/edge_list_0.tsv" \
        -s "$i" \
        -o "$intermediate/network/edge_list_${i}.tsv"
    done

echo "Permuting scores..."

    cp "$scores" "$intermediate/network_score/scores_0.tsv"

    python $PWD/src/find_permutation_bins.py \
      -gsf "$intermediate/network_score/scores_0.tsv" \
      -igf "$index" \
      -elf "$edges" \
      -ms 1000 \
      -o "$intermediate/network_score/score_bins.tsv"

    for i in `seq $num_permutations`; do
      python $PWD/src/permute_scores.py \
        -i "$intermediate/network_score/scores_0.tsv" \
        -bf "$intermediate/network_score/score_bins.tsv" \
        -s "$i" \
        -o "$intermediate/network_score/scores_${i}.tsv"
    done

################################################################################
#
#   Construct hierarchies.
#
################################################################################

echo "Constructing hierarchies..."

    for i in $(seq 0 $num_permutations); do
      python $PWD/src/construct_hierarchy.py \
        -smf "$intermediate/network/similarity_matrix.h5" \
        -igf "$index" \
        -gsf "$intermediate/network_score/scores_${i}.tsv" \
        -helf "$intermediate/network_score/hierarchy_edge_list_${i}.tsv" \
        -higf "$intermediate/network_score/hierarchy_index_gene_${i}.tsv"
    done

################################################################################
#
#   Process hierarchies.
#
################################################################################

echo "Processing hierarchies..."

# This example uses -lsb/--lower_size_bound 1 because it is a small toy example
# with 25 vertices.  Use larger value (default is 10) for larger graphs.
    python $PWD/src/process_hierarchies.py \
      -oelf "$intermediate/network_score/hierarchy_edge_list_0.tsv" \
      -oigf "$intermediate/network_score/hierarchy_index_gene_0.tsv" \
      -pelf $(for i in $(seq $num_permutations); do echo "$intermediate/network_score/hierarchy_edge_list_${i}.tsv "; done) \
      -pigf $(for i in $(seq $num_permutations); do echo "$intermediate/network_score/hierarchy_index_gene_${i}.tsv "; done) \
      -lsb 1 \
      -cf "$output_file" \
      -pl "network" "score" \
      -pf "$results/sizes_network_score.pdf"