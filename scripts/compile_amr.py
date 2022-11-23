import pandas as pd
import sys

# read in the data
orig_df = pd.read_csv(sys.argv[1]) 

# set thresholds and arguments

min_read = int(sys.argv[2])
min_read_breadth = int(sys.argv[3])
min_read_depth = int(sys.argv[4])

output_filename_root = sys.argv[5]

#
# create gene x samples matrix
#

orig_matrix = pd.pivot_table(orig_df, values=['num_reads'], index=['gene_name'], columns=['sample_name'])
orig_matrix = (orig_matrix > 0) * 1

filtered_df = orig_df.loc[(orig_df['num_reads'] >= min_read) & (orig_df['read_coverage_breadth'] >= min_read_breadth) & (orig_df['read_coverage_depth'] >= min_read_depth)]
filtered_matrix = pd.pivot_table(filtered_df, values=['num_reads'], index=['gene_name'], columns=['sample_name'])
filtered_matrix = (filtered_matrix > 0) * 1

final_matrix = orig_matrix.add(filtered_matrix, fill_value = 0) 

# in this output vvv the counts are...
# ... 0 if the gene was not present in the sample
# ... 1 if the gene was present in the sample, but didn't pass the threshold filters
# ... 2 if the gene was present in the sample and passed the threshold filters
final_matrix['num_reads'].to_csv(output_filename_root + "_gene_by_samples.csv")

# 
# create drug class x samples matrix (where value = count of genes in the drug class)
#

# function to convert data into format that separates the drug classes 
# ... and enables counting genes per class
def make_drug_class_df(input_df):
	all_rows = {}
	sum_dc = 0
	for i in input_df.index:
	    this_dict = input_df.iloc[i].to_dict()
	    drug_classes = [j.strip() for j in this_dict['drug_class'].split(';')]
	    for dc in drug_classes:
	        sum_dc += 1
	        new_dict = dict(this_dict)
	        new_dict['new_drug_class'] = dc
	        all_rows[sum_dc] = new_dict

	updated_df = pd.DataFrame.from_dict(all_rows, orient='index')
	return(updated_df)

drug_class_df = make_drug_class_df(orig_df)
drug_class_matrix = pd.pivot_table(drug_class_df, values=['gene_name'], index=['new_drug_class'], columns=['sample_name'], aggfunc='count')
drug_class_matrix_bool = (drug_class_matrix > 1000000) * 1  # to create the scaffold, but all values are set to 0

filtered_drug_class_df = drug_class_df.loc[(drug_class_df['num_reads'] >= min_read) & (drug_class_df['read_coverage_breadth'] >= min_read_breadth) & (drug_class_df['read_coverage_depth'] >= min_read_depth)]
filtered_drug_class_matrix = pd.pivot_table(filtered_drug_class_df, values=['gene_name'], index=['new_drug_class'], columns=['sample_name'], aggfunc='count')
filtered_drug_class_matrix_bool = (filtered_drug_class_matrix > 0) * 1

final_drug_class_matrix = drug_class_matrix_bool.add(filtered_drug_class_matrix, fill_value = 0)

# in this output vvv the counts correspond to the number of genes that pass the threshold and are in that drug class
final_drug_class_matrix['gene_name'].to_csv(output_filename_root + "_drug_class_by_samples.csv")