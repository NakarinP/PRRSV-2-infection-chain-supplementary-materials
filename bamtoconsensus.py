import pysam
from collections import Counter
import argparse

def bam_to_consensus_no_ref(bam_file, output_fasta, fasta_header):
    # Open the BAM file without requiring reference information
    bam = pysam.AlignmentFile(bam_file, "rb", check_sq=False)  # 'check_sq=False' avoids loading reference sequences
    
    # Dictionary to store consensus bases
    consensus = {}

    # Iterate through each pileup column (each position in the genome)
    for pileupcolumn in bam.pileup():
        # Position in 1-based coordinates
        pos = pileupcolumn.reference_pos + 1

        # Get bases at this position
        bases = [
            pileupread.alignment.query_sequence[pileupread.query_position]
            for pileupread in pileupcolumn.pileups
            if pileupread.query_position is not None
        ]

        # Count the bases
        base_counts = Counter(bases)

        # Determine the most common base (consensus)
        if base_counts:
            consensus[pos] = base_counts.most_common(1)[0][0]
        else:
            # If no base is present, use '-'
            consensus[pos] = "-"

    # Close BAM file
    bam.close()

    # Write consensus sequence to a FASTA file
    with open(output_fasta, "w") as fasta_out:
        # Write the custom header
        fasta_out.write(f">{fasta_header}\n")
        # Write the consensus sequence
        fasta_out.write("".join(consensus.get(i, "N") for i in range(1, max(consensus.keys()) + 1)))
        fasta_out.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Merge consensus sequences to non-ambiguous sequence')
    parser.add_argument('-b', '--bam', type=str, required=True, help="Read alignment in BAM format")
    parser.add_argument('-o', '--out', type=str, required=True, help="Output Fasta File Name")
    parser.add_argument('-n', '--name', type=str, required=True, help="Fasta header Name")
    args = parser.parse_args()

    bam_to_consensus_no_ref(args.bam, args.out, args.name)