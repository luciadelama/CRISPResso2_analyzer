import re

def _get_total_reads(filename):
    """Extracts the total number of reads aligned from a tab-separated file."""

    ## TOTAL READS ALIGNED ##
    infile = open(filename, 'r')

    # Read each line separately
    header = infile.readline().strip().split('\t') # List
    values = infile.readline().strip().split('\t')

    data = dict(zip(header, values))  # Dictionary
    # print(data)

    # Select the key and get the value
    total_reads = int(data["READS ALIGNED"])
    return total_reads


def _get_modif_reads(filename, total_reads):
    """Extracts the counts and percentages of modified and unmodified reads from a tab-separated file"""

    ## MODIFIED & UNMODIFIED READS ##
    infile = open(filename, 'r')

    header = infile.readline().strip().split('\t')
    values = infile.readline().strip().split('\t')

    data = dict(zip(header, values))
    # print(data)

    # Get the data
    modified_reads = int(data["Modified"])
    unmodified_reads = int(data["Unmodified"])

    wt_unmodified = round(unmodified_reads / total_reads * 100, 2)
    indel = round(modified_reads / total_reads * 100, 2)

    return {
        "modified_reads": modified_reads,
        "unmodified_reads": unmodified_reads,
        "wt_unmodified_percent": wt_unmodified,
        "indel_percent": indel,
    }


def _get_frame_reads(filename, total_reads):
    """Extracts the counts and percentages of inframe and frameshift reads from a tab-separated file"""

    ## FRAMESHIFT ANALYSIS ##
    infile = open(filename, 'r')

    for line in infile:
        line = line.strip()
        if line.startswith("In-frame mutation:"):
            inframe_reads = int(line.split(':')[1].split()[0])
        if line.startswith("Frameshift mutation:"):
            frameshift_reads = int(line.split(':')[1].split()[0])

    inframe = round(inframe_reads / total_reads * 100, 2)
    frameshift = round(frameshift_reads / total_reads * 100, 2)

    return {
        "inframe_reads": inframe_reads,
        "frameshift_reads": frameshift_reads,
        "inframe_percent": inframe,
        "frameshift_percent": frameshift,
    }

# def _get_input_run(filename):
    """Extracts --amplicon_seq and --coding_seq values from a CRISPResso command line in a log file."""

    cmd_line = None
    amplicon_seq = None
    coding_seq = None

    with open(filename, 'r') as infile:
        for line in infile:
            if line.startswith("[Command used]:"):
                cmd_line = next(infile).strip()
                parts = cmd_line.split()

                for i, part in enumerate(parts):
                    if part in ("--amplicon_seq", "-a") and i + 1 < len(parts):
                        amplicon_seq = parts[i + 1]
                    elif part in ("--coding_seq", "-c") and i + 1 < len(parts):
                        coding_seq = parts[i + 1]

                if amplicon_seq or coding_seq:
                    break
    
    return {
        "amplicon_seq": amplicon_seq,
        "coding_seq": coding_seq
    }

def _get_mut_wt_reads(filename, wt_seq, mut_seq, amplicon, total_reads):
    """Reads a tab-separated file and calculates mutated and wild-type read counts and percentages."""

    wt_reads = None
    mut_reads = None

    with open(filename, 'r') as infile:
        for line in infile:
            parts = line.strip().split('\t')
            aligned_seq = parts[0]
            ref_seq = parts[1]

            if wt_reads is None and re.search(aligned_seq, wt_seq) and re.search(ref_seq, amplicon):
                wt_reads = int(parts[6])
            
            elif mut_reads is None and re.search(aligned_seq, mut_seq) and re.search(ref_seq, amplicon):
                mut_reads = int(parts[6])

            if wt_reads is not None and mut_reads is not None:
                break

    mut_percent = round(mut_reads/total_reads * 100, 2)
    wt_percent = round(wt_reads/total_reads * 100, 2)
    mut_wt_percent = round(mut_reads/wt_reads * 100, 2)

    return {
        "mut_reads": mut_reads,
        "wt_reads": wt_reads,
        "mut_percent": mut_percent,
        "wt_percent": wt_percent,
        "mut_wt_percent": mut_wt_percent,
    }

# def _get_sensitivity(num_value, denom_value):
    sensitivity_value = num_value/denom_value
    return sensitivity_value