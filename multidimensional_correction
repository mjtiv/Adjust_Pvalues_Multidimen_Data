#!/usr/bin/env python3.6

import numpy as np
import sys


"""

PROGRAM: Multidimensional p-value Correction

DESCRIPTION: Program takes in a multi-dimensional dataset (see example file)
and begins correcting p-values across the entire dataset, taking into account
the multidimensial directionality of the data, so individual sample p-values for
a specific variant can be examined for significance. 

Program Based on the Paper:

Controlling False Discoveries in Multidimensional Directional
Decisions, with Applications to Gene Expression Data
on Ordered Categories
Wenge Guo,1 Sanat K. Sarkar,2 and Shyamal D. Peddada3


"""


def convert_pvalues_to_floats(samples_pvalues):

    """
    Converts a list of string values to floats and removes
    all NaNs

    :param samples_pvalues: list of values

    :return list_of_floats: list of floats all NaNs removed

    """

    # Create an empty list to store values
    list_of_floats = []

    # Start looping over the data
    for value in samples_pvalues:

        # Convert the value to a float
        value = float(value)

        # Determine if the value contains an NaN
        # Skip if it does
        if np.isnan(value)== True:
            continue
        
        # Add all passing values to the list
        else:
            list_of_floats.append(value)

    return (list_of_floats)
    


def create_pvalues_fdr_results_dict(pvalues_list):

    '''

    Creates a dictionary inside a list to allow corrected
    FDR pvalues to be stored and also allow for values to changed
    and updated in the inside dictionary

    : parameters pvalues_list: list of input pvalues

    : return pvalues_fdr_results_list: list dictionary where the list
        value is the original p-value

    '''

    # Create list dictionary combo to store results
    pvalues_list_dict = []

    # Start index counter
    x = 0

    # Start looping over pvalues
    for pvalue in pvalues_list:

        # For each entry add a value as an outside list with a dictionary inside
        # Allow sorting of the outside value
        pvalues_list_dict.append([pvalue, {'index_in_list': x,
                                               'original_pvalue': pvalue,
                                               'corrected_pvalue': 'nan'}])
        # Add to the index counter
        x += 1

    return (pvalues_list_dict)



def reorder_sorted_pvalues_list_dict(sorted_pvalues_list_dict):

    """

    Takes in the final list-dictionary with corrected p-values
    and replaces pvalue "list" vlue with the original list index value
    so the list-dictionary can be re-sorted to match the original list
    order.

    : parameters sorted_pvalues_list_dict: list-dictionary with all entries
        where the list value is the original p-value

    : return re_sorted_pvalues_list_dict: list-dictionary where the
        original the list value is the original index value and has been
        re-sorted by that value

    """

    for entry in sorted_pvalues_list_dict:

        index_value = entry[1]['index_in_list']

        entry[0] = index_value

    # Sort the list (outside value using original index value)
    re_sorted_pvalues_list_dict = sorted(sorted_pvalues_list_dict)

    return (re_sorted_pvalues_list_dict)



def fdr_correction(pvalues_list):

    """

    Function performs FDR correction of the p-values
    using Benjamini-Hochberg (1995) which sorts the list of pvalues
    and then determines the p-value correction based on the rank and following
    equation (p-value x NumbTest / p-value_Rank) and also adjusts for lower p-value 
    in succeeding adjusted p-value in list if it occurs.

    :param values: list of pvalues to correct

    :return p_value_dict: dictionary of corrected pvalues
    
    """
    
    # Create a list of pvalues with a dictionary inside each entry
    pvalues_list_dict = create_pvalues_fdr_results_dict(pvalues_list)

    # Sort the list (outside value used raw pvalue)
    sorted_pvalues_list_dict = sorted(pvalues_list_dict)
    
    #Position Movement Counter
    i=1

    total_pvalues = len(pvalues_list_dict)

    # Start looping over list dictionary 
    for entry in sorted_pvalues_list_dict:

        if i < len(pvalues_list_dict):

            # Get corrected pvalue for position i
            adj_pvalue = round((entry[0] * total_pvalues / i), 8)
            
            # Get next pvalue in sorted list
            next_entry = sorted_pvalues_list_dict[i][0]

            # Get corrected pvalue for next entry
            next_entry_adj_pvalue = round((next_entry * total_pvalues / (i + 1)), 8)

            # Increment i value
            i+=1

            # Get the lowest possible value of the two options
            fdr_corrected_pvalue = min(adj_pvalue, next_entry_adj_pvalue)

            # Update the dictionary value with corrected pvalue
            entry[1]['corrected_pvalue'] = fdr_corrected_pvalue

        # last value in list
        else:
 
            # Get corrected pvalue for position i
            # Important: No correction occurs for position because last value in sorted list
            fdr_corrected_pvalue = round((entry[0] * total_pvalues / i), 8)

            # Update the dictionary value with corrected pvalue
            entry[1]['corrected_pvalue'] = fdr_corrected_pvalue

    # Replace first value in list with index to re-order values to match orginal data
    re_sorted_pvalues_list_dict = reorder_sorted_pvalues_list_dict(sorted_pvalues_list_dict)

    return(re_sorted_pvalues_list_dict)


def determine_passing_pvalues(lowest_pvalue_list,
                              fdr_corrected_dict):

    """

    Loops over the dictionary of pvalues using the original
    uncorrected pvalue list to control order and determines
    all FDR corrected pvalues that pass the cutoff

    : param lowest_pvalue_list: list of lowest pvalues
    : fdr_corrected_dict: dictionary of all corrected pvalues

    : return count_passing_pvalues: count of all pvalues that pass the cutoff

    """

    # List to store all passing pvalues
    passing_pvalues = []

    # Inex in List Counter for Movement
    index_in_list_counter = 0

    # Looping over the dictionary
    for value in lowest_pvalue_list:

        # Get the dictionary results part of the list-dictionary combo
        dictionary_results =(fdr_corrected_dict[index_in_list_counter][1])

        # Add to index_in_list_counter now that values have been retrieved
        index_in_list_counter += 1

        # Get the original pvalue for double-checking it matches
        original_pvalue = dictionary_results['original_pvalue']

        # Get the adjusted pvalue
        fdr_pvalue = dictionary_results['corrected_pvalue']

        # Safety to kill program if lists being sorted over index order is wrong
        if original_pvalue != value:
            print ("P-value Lists Order is not correct, double check code")
            print ("Killing program")
            sys.exit(1)

        # Determine if pvalue is below the 0.05 cutoff
        if fdr_pvalue < 0.05:

            # Add to the list if it passes
            passing_pvalues.append(fdr_pvalue)

        else:
            continue

    # Get the count of passing values
    count_passing_pvalues = len(passing_pvalues)

    return (count_passing_pvalues)
        

def determine_passing_FDR_pvalues(input_file_name):

    """
    Opens the inputted file and reads through the data extracting
    the lowest pvalue for each variant and Bonferroni corrects. The final
    pvalue is then FDR corrected and then the final list of FDR corrected
    pvalues examined to see how many pass 0.05 cutoff.

    : parameter input_file_name: name of input file to open and analyze

    : return passing_count: number of variants that pass overall correction

    """

    # Open the file
    input_file = open(input_file_name, 'r')

    # Line movement counter
    line_counter = 0

    # Lowest p-value list
    lowest_pvalue_list = []

    # Total variants count
    variant_counter = 0

    # Read through the lines of the file
    for line in input_file:

        # Skip over the header line
        if line_counter == 0:
            line_counter += 1
            continue

        # Start parsing variant results
        else:

            # Start counting variants
            variant_counter +=1

            # Clean up line
            line = line.rstrip("\n")

            # Split the data
            data = line.split("\t")

            # Seperate out the parts of the line
            variant_id = data[0]
            samples_pvalues = data[1:]

            # Convert values to floats and also remove NaNs
            samples_pvalues = convert_pvalues_to_floats(samples_pvalues)

            # Find lowest p-value
            min_pvalue = np.amin(samples_pvalues)

            # Length of list
            length_of_list = len(samples_pvalues)

            # Calculate the bonferroni pooled p-value
            bonf_pooled_pvalue = min_pvalue * length_of_list
            
            # Add value to list for FDR correction
            lowest_pvalue_list.append(bonf_pooled_pvalue)

    # Get FDR corrected p-values
    fdr_corrected_dict = fdr_correction(lowest_pvalue_list)

    # Determine the number of passing p-values after FDR correction
    passing_count = determine_passing_pvalues(lowest_pvalue_list,
                                              fdr_corrected_dict)

    # Close the input file
    input_file.close()

    return(passing_count, variant_counter)

          
def determine_passing_variant_pvalues(input_file_name, passing_variant_count, total_variants_analyzed):

    """

    Opens up the input_file and starts parsing through all the results determining if
    the p-values pass the new corrected p-value based on the entire dataset (multi-dimensional correction).
    Outputs a new file with whether the prior p-value passes the new multi-dimensional corrected cutoff value

    : parameter input_file_name: name of the file being analyzed
    : passing_variant_count: number of variants that pass the prior correction
    : total_variants_analyzed: number of variants analyzed in the original dataset

    : return NONE:

    """

    # Open the file
    input_file = open(input_file_name, 'r')

    # Output file
    output_file = open("corrected_test_data.txt", "w")

    # Line movement counter
    line_counter = 0

    # Lowest p-value list
    lowest_pvalue_list = []

    # Read through the lines of the file
    for line in input_file:

        # Skip over the header line
        if line_counter == 0:
            line_counter += 1

            # Write header to new file
            output_file.write(line)
            continue

        # Start parsing variant results
        else:
            line = line.rstrip("\n")
            data = line.split("\t")

            # Seperate out the parts of the line
            variant_id = data[0]
            samples_pvalues = data[1:]

            # Print variant name to file
            output_file.write(variant_id + "\t")

            # Convert values to floats
            number_of_tests = len(convert_pvalues_to_floats(samples_pvalues))

            # Calculate new p-value cutoff to compare to original p-value (round to 9 decimal places)
            # Note: The rounding prevents issues with how python stores trailing zeroes
            corrected_pvalue_threshold = round((passing_variant_count * 0.05) / (number_of_tests * total_variants_analyzed), 9)

            # Start looping over original values and correcting
            sample_counter = 1
            for value in samples_pvalues:

                # Identify last sample in list to print newline and not tab
                if sample_counter == len(samples_pvalues):

                    # Convert the value to a float
                    value = float(value)

                    # Determine if the value contains an NaN
                    # Skip if it does
                    if np.isnan(value)== True:
                        output_file.write(str(value) + "\n")

                    # Non NaN values start analyzing
                    else:
                        # Analyze value to see if passes new threshold
                        if value < corrected_pvalue_threshold:
                            output_file.write('pass|' + str(value) + "\n")

                        # If value fails the new cutoff threshold
                        else:
                            output_file.write('fail|' + str(value) + "\n")

                else:
                    # Add to sample counter
                    sample_counter +=1
                    
                    # Convert the value to a float
                    value = float(value)

                    # Determine if the value contains an NaN
                    # Skip if it does
                    if np.isnan(value)== True:
                        output_file.write(str(value) + "\t")

                    # Non NaN values start analyzing
                    else:

                        # Analyze value to see if passes new threshold
                        if value < corrected_pvalue_threshold:
                            output_file.write('pass|' + str(value) + "\t")

                        # If value fails the new cutoff threshold
                        else:
                            output_file.write('fail|' + str(value) + "\t")

    # Close the file
    output_file.close()
      
    return()


def main():

    print ("Starting Multi_Dimensional Pvalue Adjustment Program")

    input_file_name = 'test_data.txt'

    # Determine number of pass FDR pvalues after Bonferroni pooling and FDR
    passing_FDR_pvalue_results = determine_passing_FDR_pvalues(input_file_name)

    # Get actual results from analysis
    passing_variant_count = passing_FDR_pvalue_results[0]
    total_variants_analyzed = passing_FDR_pvalue_results[1]

    # Determine if pvalues for each variant pass the adjustment
    determine_passing_variant_pvalues(input_file_name, passing_variant_count, total_variants_analyzed)

    print ("Program Ran Correctly")

if __name__ == "__main__":
    
    main()








