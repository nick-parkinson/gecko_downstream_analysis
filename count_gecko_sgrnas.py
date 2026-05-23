# -*- coding: utf-8 -*-
"""
Count sgrna reads in fastq files from CRISPR screen.
Outputs in correct format for MAGeCK

Assumes all sgRNAs are 20nt in length.

Requires fastq files, and list of guide seqs and targets

Guide seq / target file format:
    tab-delimited file, column 1 = gene, col 2 = seq, no header.

Will output single file of raw counts for all pools/conditions,
+ summary file of library stats for each pool.

Created on Thu Oct 24 11:20:03 2019

@author: nparkins
"""

import pandas as pd
import numpy as np
pd.set_option('mode.chained_assignment', None)

# pool_id: pool fastq file path. Note, 4 files per pool.
POOLS = {
         'pool1_name': '/path/to/pool1_fastqfile.fastq',
         'pool2_name': '/path/to/pool2_fastqfile.fastq'

        }

GUIDELIST = "/path/to/sgrna_seq_target_file.txt"

OUTPUTFILE = "/path/to/outputfile.txt"

# Output details of 20-nt guides not matching library, for troubleshooting
OUTPUT_NONMATCHING = False  

# Vector-specific protospacer flank sequences:
# Sequence immediately preceeding protospacer:
KEY1 = "CGAAACACCG"  

# expected sequence following protospacer
# (e.g. GTTTT for LentiCRISPRv2, GTTT for pkLV2)
KEY2 = "GTTT"  

# Columns for guide list:
SEQ_COLUMN = 1
GENE_COLUMN = 0

MIN_QUAL = ','  # lowest acceptable fastq quality score

# ------------------- FUNCTIONS -------------------------------------------
def get_guides(file):
    """
    Makes dictionary of guide seq: target
    """
    with open(file) as f:
        lines = [i.split("\t") for i in f.readlines()]
        # check for header:
        if len(lines[0][SEQ_COLUMN]) < 20:
            guide_dict = {line[SEQ_COLUMN].strip("\n").upper():
                      line[GENE_COLUMN].strip("\n") for line in lines[1:]}
        else:
            guide_dict = {line[SEQ_COLUMN].strip("\n").upper():
                      line[GENE_COLUMN].strip("\n") for line in lines}        
        return guide_dict


def get_sgrna(read, qual, keylength, key1, key2, minq):
    """
    Finds guide seq from read,
    if quality is sufficient.
    Ignore read if ANY base is below threshold score
    """
    notfound=0
    wronglength=0
    lowqual = 0
    guide=False
    findkey = read.find(key1)
    if findkey == -1:
        notfound = 1
    else:
        if not read[findkey+20+keylength:].startswith(key2):
            wronglength = 1
            if OUTPUT_NONMATCHING:
                # get 20mer for error report
                guide = read[findkey+keylength:findkey+keylength+20]
        elif any([ord(i) < minq for i in qual[findkey+keylength:findkey+keylength+20]]):
            lowqual = 1
        else:
            guide = read[findkey+keylength:findkey+keylength+20]
            
    return guide, notfound, wronglength, lowqual


def count_sgrnas(read_list, guide_list):
    """
    Iterates through fastq reads, looking for 20-nt sgRNA sequences between
    'key' sequences. Logs counts for guides and non-matching sequences.
    In this version, will amalgamate different lane files for the same sample
    """
    guides = {}
    notfound_count = 0  # reads where key not found
    wronglength_count = 0  # keys not separated by 20-nt
    lowqual_count = 0  # read quality below threshold for 1 or more bases of guide seq
    keylength = len(KEY1)
    minq = ord(MIN_QUAL)  # convert character to ascii int
    readcounter = 0
    # Up to 4 lanes...
    if 'L001' in read_list:
        filelist = [read_list.replace("L001", "L00{}".format(l)) for l in [1,2,3,4]]
    else:
        filelist = [read_list]
    for file in filelist:
        try:
            with open(file) as f:
                counter = 0
                for line in f:
                    if counter == 1:
                        thisread = line.strip("\n")
                        readcounter += 1
                    if counter == 3:
                        thisqual = line.strip("\n")
                        guide, notfound, wronglength, lowqual = get_sgrna(thisread,
                            thisqual, keylength, KEY1, KEY2, minq)
                        if not guide:
                            notfound_count += notfound
                            wronglength_count += wronglength
                            lowqual_count += lowqual                      
                        else:
                            if OUTPUT_NONMATCHING and (wronglength == 1):
                                wronglength_count += wronglength
                            if guide not in guides:
                                guides[guide] = 1
                            else:
                                guides[guide] += 1
                        counter = 0
                    else:
                        counter += 1
        except FileNotFoundError:
            continue
    # fill in undetected guides from expected guide list
    for guide in guide_list:
        if guide not in guides:
            guides[guide] = 0
    print(readcounter)
    return guides, notfound_count, wronglength_count, lowqual_count, readcounter


def filter_read_errors(guide_counts, guide_list):
    """
    Splits guide count dictionary into matching and non-matching guides
    """
    true_dict = {}
    error_dict = {}
    for key, val in guide_counts.items():
        if key in guide_list:
            true_dict[key] = val
        else:
            error_dict[key] = val
    return true_dict, error_dict


def calculate_stats(true, errors, notfound, wrong_length, low_qual, guide_list, readcounter):
    """
    Calculate library statistics
    """
    stats = {}
    stats['total_reads'] = readcounter
    stats['perfect_matches'] = sum(true)
    stats['key_not_found'] = notfound
    stats['wrong_length'] = wrong_length
    if len(errors) == 0:
        stats['nonmatching'] = 0
    else:
        if OUTPUT_NONMATCHING:
            # wrong length will be included in errors
            # nonmatching gives 20mers only
            stats['nonmatching'] = sum(errors) - wrong_length
        else:
            stats['nonmatching'] = sum(errors)
    stats['low_quality'] = low_qual
    stats['all_incorrect'] = stats['nonmatching'] + wrong_length
    #stats['total_reads'] = stats['perfect_matches'] + stats['all_incorrect'] + notfound + low_qual
    stats['perc_key_not_found'] = 100 * stats['key_not_found']/stats['total_reads']
    stats['perc_wronglength'] = 100 * stats['wrong_length']/stats['total_reads']
    stats['perc_nonmatching'] = 100 * stats['nonmatching']/stats['total_reads']
    stats['perfect_perc'] = 100 * stats['perfect_matches']/stats['total_reads']
    stats['perfect_perc_hqual'] = 100 * stats['perfect_matches']/(stats['total_reads']-low_qual)
    stats['incorrect_perc'] = 100 * stats['all_incorrect']/stats['total_reads']
    stats['lq_perc'] = 100 * stats['low_quality']/stats['total_reads']
    stats['total_guides'] = len(guide_list)
    stats['guides_not_found'] = list(true).count(0)
    stats['guides_found'] = stats['total_guides'] - stats['guides_not_found']
    stats['perc_guides_found'] = 100 * stats['guides_found'] / stats['total_guides']
    stats['upper'] = np.percentile(true, 90)
    stats['lower'] = np.percentile(true, 10)
    if stats['upper'] != 0 and stats['lower'] != 0:
        stats['skew_ratio'] = stats['upper']/stats['lower']
    else:
        stats['skew_ratio'] = 'Not enough perfect matches to determine skew ratio'

    return stats


def write_summary(name, stats):
    """
    Write summary text file with library stats for each pool
    """
    with open("pool_{}_summary.txt".format(name), 'w') as o:
        o.write("Total reads: {}\n".format(stats['total_reads']))
        o.write("Key not found: {} reads ({}%)\n".format(
                stats['key_not_found'], stats['perc_key_not_found']))
        o.write("Wrong length spacer: {} reads ({}%)\n".format(stats['wrong_length'], stats['perc_wronglength']))
        o.write("Unmatched 20mers: {} reads ({}%)\n".format(stats['nonmatching'], stats['perc_nonmatching']))
        o.write("Incorrect guide sequence: {} reads ({}%)\n".format(
                stats['all_incorrect'], stats['incorrect_perc']))
        o.write("Low read quality, read ignored: {} ({}%)\n".format(
                 stats['low_quality'], stats['lq_perc']))
        o.write("Reads perfectly matching guides: {}, ({}% of total, {}% of high quality reads)\n".format(
                stats['perfect_matches'], stats['perfect_perc'],
                stats['perfect_perc_hqual']))
        o.write("Expected guides found: {}/{} ({}%)\n".format(
                stats['guides_found'], stats['total_guides'],
                stats['perc_guides_found']))
        o.write("90th percentile count: {}\n".format(
                stats['upper']))
        o.write("10th percentile count: {}\n".format(
                stats['lower']))
        o.write("Skew ratio: {}".format(stats['skew_ratio']))


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Get dictionary of expected guides and their targets
    guide_dict = get_guides(GUIDELIST)
    guide_list = set(guide_dict)

    # Create dataframe for count results, sequence as index
    df = pd.DataFrame.from_dict(guide_dict, orient='index', columns=['target'])

    # For quality control / troubleshooting, also store details of 20-nt
    # guide sequences that do not perfectly match a guide in the library
    df_errors = pd.DataFrame()

    for pool, file in POOLS.items():
        # count guides + non-matching reads
        guide_counts, notfound, wrong_length, low_quality, readcounter = count_sgrnas(file, guide_list)

        # Identify reads matching known guides or not, add to appropriate df
        true_dict, error_dict = filter_read_errors(guide_counts, guide_list)
        if len(error_dict.keys()) == 0:
            print("No errors found for pool {}".format(pool))
        df[pool] = pd.Series(true_dict)
        error_temp = pd.DataFrame.from_dict(error_dict, orient='index')
        error_temp.columns = [pool]
        df_errors = pd.concat([df_errors, error_temp], axis=1,
                              copy=False, sort=True).fillna(0)

        # Calculate and report library stats for this pool
        stats_dict = calculate_stats(df[pool].values, df_errors[pool].values,
                                     notfound, wrong_length, low_quality, guide_list, readcounter)
        write_summary(pool, stats_dict)

    # Output counts of true guides and errors if wanted
    df.to_csv(OUTPUTFILE, sep="\t", index_label='seq')
    if OUTPUT_NONMATCHING:
        df_errors.to_csv("errors.{}".format(OUTPUTFILE), sep="\t")
