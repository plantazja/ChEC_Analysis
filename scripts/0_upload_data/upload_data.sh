#!/bin/bash

# Parse command line arguments
usage () {
    echo "Usage: bash upload_data.sh -i /s/donczew-lab/RDHTSP* -o /hpc-prj/donczew/Anastasia/ChEC_analysis/Data"
}

while getopts "i:o:h" opt;
do
    case $opt in
        i) 
            in_dir=${OPTARG}
            ;;
        o) 
            out_dir=${OPTARG}
            ;;
        h) 
            usage 
            exit 0
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "$in_dir" ] || [ -z "$out_dir" ]; then
    echo "Error: Both -i and -o arguments are required"
    usage
    exit 1
fi

cd "$(dirname "$0")"

# Transfer bam files
bam_in_dir="${in_dir}/ChEC-Seq_Analysis-*/results/alignment/scer"

mkdir -p "$out_dir"
bam_out_dir="${out_dir}/bam"
mkdir -p "$bam_out_dir"

for file in $bam_in_dir/*.bam
do     
    filename=$(basename "$file")
    if [ -f "${bam_out_dir}/${filename}" ]; then
        echo "BAM file already exists in $bam_out_dir, skipping.."
    else
        echo "Copying $filename BAM file to $bam_out_dir.."
        cp "$file" "$bam_out_dir/"
    fi
done

# Transfer bigwig files
bwig_out_dir="${out_dir}/bigwig"
mkdir -p "$bwig_out_dir"

cp -r "${in_dir}"/ChEC-Seq_Analysis*/results/bigwig/cpm "$bwig_out_dir/" 
cp -r "${in_dir}"/ChEC-Seq_Analysis*/results/bigwig/cpm_mean "$bwig_out_dir/" 
cp -r "${in_dir}"/ChEC-Seq_Analysis*/results/bigwig/spikein "$bwig_out_dir/" 
cp -r "${in_dir}"/ChEC-Seq_Analysis*/results/bigwig/spikein_mean "$bwig_out_dir/"

# Convert bigwig to wig
BIGWIG_TOOL="./bigWigToWig" 
wig_out_dir="${out_dir}/wig"
mkdir -p "$wig_out_dir"

# Convert each .bw file
for bw_subdir in "$bwig_out_dir"/*/ ; do
    # Only process if it's a directory
    if [ -d "$bw_subdir" ]; then
        subdir_name=$(basename "$bw_subdir")
        
        # Create corresponding subdirectory in wig_out_dir
        mkdir -p "$wig_out_dir/$subdir_name"
        
        # Convert each .bw file in the subdirectory
        for bw_file in "$bw_subdir"/*.bw; do
            if [ -f "$bw_file" ]; then
                filename=$(basename "$bw_file" .bw)
                wig_file="$wig_out_dir/$subdir_name/$filename.wig"
                
                if [ -f "$wig_file" ]; then
                    echo ".wig file already exists in $wig_out_dir/$subdir_name, skipping.."
                else
                    echo "Converting $filename.bw to wig in $subdir_name..."
                    "$BIGWIG_TOOL" "$bw_file" "$wig_file"
                    
                    if [ $? -eq 0 ]; then
                        echo "Done"
                    else
                        echo "Failed to convert $bw_file"
                    fi
                fi
            fi
        done
    fi
done

echo "Upload complete!"