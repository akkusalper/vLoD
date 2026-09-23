import pandas as pd
import sys

def merge_detectability_into_vcf(vcf_path, detectability_path, output_path):
    # Read the detectability status file
    xls_df = pd.read_csv(detectability_path, sep='\t')
    detectability_data_dict = {}
    for index, row in xls_df.iterrows():
        detectability_data_dict[row['VCF_ID']] = {
            'condition': 'Yes' if row['Detectability_Condition'] == 'Detectable' else 'No',
            'score': row['Detectability_Score']
        }

    # Process the VCF file
    updated_vcf_lines = []
    with open(vcf_path, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("##INFO"):
                updated_vcf_lines.append(line)
                if "TLOD" in line:
                    updated_vcf_lines.append('##INFO=<ID=DET,Number=A,Type=String,Description="Detectability status per ALT allele (Yes if detectable, No if non-detectable)">\n')
                    updated_vcf_lines.append('##INFO=<ID=DETS,Number=A,Type=Float,Description="Detectability Score per ALT allele">\n')
                continue
            if line.startswith("##") or line.startswith("#"):
                updated_vcf_lines.append(line)
                continue
            columns = line.strip().split("\t")
            chrom, pos, ref, alts = columns[0], columns[1], columns[3], columns[4]

            # The detectability scorer emits one row per decomposed ALT allele
            # (e.g. chr12_25245351_C_A and chr12_25245351_C_G), so a
            # multiallelic VCF line's raw, comma-joined ALT field ("A,G")
            # never matches a dict key on its own. Look up each allele
            # separately and recombine, per VCF's Number=A convention.
            statuses, scores = [], []
            for alt in alts.split(","):
                entry = detectability_data_dict.get(f"{chrom}_{pos}_{ref}_{alt}")
                statuses.append(entry['condition'] if entry else ".")
                scores.append(str(entry['score']) if entry else ".")

            if any(status != "." for status in statuses):
                columns[7] = columns[7] + f";DET={','.join(statuses)};DETS={','.join(scores)}"
            updated_vcf_lines.append("\t".join(columns) + "\n")

    # Save the updated VCF
    with open(output_path, "w") as out_vcf:
        out_vcf.writelines(updated_vcf_lines)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_detectability.py <path_to_vcf> <path_to_detectability_status> <output_path>")
        sys.exit(1)

    vcf_path = sys.argv[1]
    detectability_path = sys.argv[2]
    output_path = sys.argv[3]
    
    merge_detectability_into_vcf(vcf_path, detectability_path, output_path)