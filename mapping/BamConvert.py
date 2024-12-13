import mapping.TsUtils as tsUtils

#
import pysam
from Bio import SeqIO


def _convertPos(pos, tsInfo):
    relpos = 0
    for s, e in tsInfo.exons:
        if pos < relpos + (e - s):
            return s + (pos - relpos)
        else:
            relpos = relpos + (e - s)

    return 0


def convertPos(r_st, r_en, tsInfo):
    if tsInfo.strand == "-" or tsInfo.strand == "-1":
        tlen = 0
        for s, e in tsInfo.exons:
            tlen = tlen + (e - s)

        r_stm = tlen - r_en
        r_edm = tlen - r_st

        r_st = _convertPos(r_stm, tsInfo)
        r_en = _convertPos(r_edm, tsInfo)

    else:

        r_en = _convertPos(r_en, tsInfo)

    r_st += 1
    return r_st, r_en


def checkN_Start(r_st_g, cigar_g):
    a = pysam.AlignedSegment()
    a.cigarstring = cigar_g
    cigarlist = []
    cnt = 0
    isCond = False
    for v in a.cigar:

        cigaroprator, cigarlen = v
        if ((cigaroprator == 3) and (cnt == 0)):  # N
            r_st_g = r_st_g - cigarlen
            isCond = True
        else:
            cigarlist.append(v)

        cnt += 1

    if isCond:
        return r_st_g, toCigarStr(cigarlist)
    else:
        return r_st_g, cigar_g  # as it is


import mappy as mp


def convert(ref, start, tsData, cigar, seq):
    # print(ref,start,cigar,seq)

    tsInfo = tsData[ref]
    chrom_g = tsInfo.chrom
    strand = True
    containHardClip = False

    cigarlist = []
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refconsume = 0
    for v in a.cigar:

        cigaroprator, cigarlen = v

        if cigaroprator == 5:
            containHardClip = True

        if cigaroprator == 0 or cigaroprator == 2:  # match or del
            refconsume = refconsume + cigarlen
        if cigaroprator == 2 and cigarlen >= 200:
            v = 3, cigarlen

        cigarlist.append(v)

    if tsInfo.strand == "-" or tsInfo.strand == "-1":
        ##rev strand
        strand = False

        end = start + refconsume
        tlen = 0
        for s, e in tsInfo.exons:
            tlen = tlen + (e - s)

        r_stm = tlen - end
        r_st_g = _convertPos(r_stm, tsInfo)
        seq_g = mp.revcomp(seq)
        cigarlist.reverse()
        cigar_g = cigarlist


    else:
        r_st_g = _convertPos(start, tsInfo)
        seq_g = seq
        cigar_g = cigarlist

    # insertN
    debug = False
    cigar_g = insertN(cigar_g, r_st_g, tsInfo, debug)
    cigar_g = mergeN_D(cigar_g)

    # check N start
    if "N" in cigar_g:
        r_st_g, cigar_g = checkN_Start(r_st_g, cigar_g)

    return chrom_g, r_st_g, strand, cigar_g, seq_g, containHardClip


def mergeN_D(cigar):
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    modified_cigar = []
    prev_op, prev_length = None, 0

    DorN = [2, 3]
    consective = False
    for op, length in a.cigar:
        # Merge consecutive 2 (D) and 3 (N) operations
        if (prev_op in DorN and op in DorN):
            prev_length += length
            consective = True
        else:
            if prev_op is not None:
                modified_cigar.append((prev_op, prev_length))
            prev_op, prev_length = op, length

    # Add the last element
    if prev_op is not None:
        modified_cigar.append((prev_op, prev_length))

    if consective:
        return toCigarStr(modified_cigar)
    else:
        return cigar  # unchange


def insertN(cigar, r_st_g, tsInfo, debug):
    start = r_st_g
    # findout N position as list
    nlist = []
    exons = tsInfo.exons
    exonLen = len(exons)
    nlensum = 0
    for n in range(1, exonLen):

        s0, e0 = exons[n - 1]
        s1, e1 = exons[n]
        pos = e0 - start - nlensum

        nlen = s1 - e0
        # print("s0,e0,s1,e1", s0, e0, s1, e1,nlen)
        if pos > 0:
            nlist.append((pos, nlen))
            nlensum = nlensum + nlen

    if len(nlist) > 0:
        return _insertN(cigar, nlist, debug)
    else:
        return toCigarStr(cigar)


def _insertN(cigar, nlist, debug):
    # Insert N
    cigarlist = []
    # print("cigar b4",cigar)
    for cigaroprator, cigarlen in cigar:

        if cigaroprator != 3:
            cigarlist.append((cigaroprator, cigarlen))

    for pos, nlen in nlist:
        if debug:
            print("insert N in", pos, nlen)
            print(toCigarStr(cigarlist))
        cigarlist = insertN_each(pos, nlen, cigarlist)

    return toCigarStr(cigarlist)


def toCigarStr(cigarlist):
    cgstr = ""
    for cigaroprator, cigarlen in cigarlist:

        if cigaroprator == 3:  # N
            cgstr = cgstr + str(cigarlen) + "N"
        elif cigaroprator == 4:  # S
            cgstr = cgstr + str(cigarlen) + "S"
        elif cigaroprator == 5:  # H
            cgstr = cgstr + str(cigarlen) + "H"
        elif cigaroprator == 6:  # S
            cgstr = cgstr + str(cigarlen) + "P"
        elif cigaroprator == 0:  # match M
            cgstr = cgstr + str(cigarlen) + "M"
        elif cigaroprator == 2:  # Del
            cgstr = cgstr + str(cigarlen) + "D"
        elif cigaroprator == 1:  # Ins
            cgstr = cgstr + str(cigarlen) + "I"
    return cgstr


def insertN_each(pos, nlen, cigarlist):
    newlist = []
    refpos = 0
    cigaropratorN = 3
    beforeadd = True
    for cigaroprator, cigarlen in cigarlist:

        if cigaroprator == 3:  # N

            newlist.append((cigaroprator, cigarlen))

        elif cigaroprator == 0 or cigaroprator == 2:  # match or Del

            if (refpos <= pos) and (pos <= refpos + cigarlen) and beforeadd:

                left = pos - refpos
                if left > 0:
                    newlist.append((cigaroprator, left))
                newlist.append((cigaropratorN, nlen))
                # print("insertn",pos, nlen)

                right = cigarlen - left
                if right > 0:
                    newlist.append((cigaroprator, right))
                beforeadd = False
            else:
                newlist.append((cigaroprator, cigarlen))

            refpos = refpos + cigarlen



        else:  # Ins or Softclip or

            newlist.append((cigaroprator, cigarlen))

    return newlist


def calculate_read_length(cigar):
    a = pysam.AlignedSegment()
    a.cigarstring = cigar

    read_length = 0
    for op, length in a.cigar:
        if op in [0, 1, 4]:  # 0: M, 1: I, 4: S
            read_length += length
    return read_length


# def reverse_numbers_in_section(section):
#     """ Reverse the order of numbers in a section of the string and return the count of numbers. """
#     parts = section.split(',')
#     numbers = [part for part in parts if part.isdigit()]
#     numbers_reversed = numbers[::-1]
#     num_iter = iter(numbers_reversed)
#
#     reversed_section = ','.join([next(num_iter) if part.isdigit() else part for part in parts])
#     number_count = len(numbers)
#
#     return reversed_section, number_count


# def reverse_subsections_of_numbers_return_list(numbers, counts):
#     """ Reverse subsections of the number sequence based on the counts provided, and return as a list of numbers. """
#     start_index = 0
#     reversed_subsections = []
#
#     for count in counts:
#         subsection = numbers[start_index:start_index + count]
#         reversed_subsection = subsection[::-1]
#         reversed_subsections.extend(reversed_subsection)
#         start_index += count
#
#     # Convert the elements to integers
#     return [int(num) for num in reversed_subsections]


# def reverse_numbers_in_string(s):
#     """ Reverse the order of numbers in each section of the string separated by semicolons and return counts. """
#     sections = s.split(';')
#     reversed_sections_info = [reverse_numbers_in_section(section) for section in sections if section]
#
#     reversed_sections = [info[0] for info in reversed_sections_info]
#     number_counts = [info[1] for info in reversed_sections_info]
#
#     return ';'.join(reversed_sections), number_counts

from array import array


class Converter:

    def __init__(self, reflist, convertMatrix):

        self.reflist = reflist
        self.tsData = tsUtils.getTsDict(convertMatrix, 1)

    def convert_read(self, read):

        cs = read.cigarstring
        ref = read.reference_name

        start = read.reference_start
        seq = read.seq
        #
        chrom_g, r_st_g, strand_g, cigar_g, seq_g, containHardClip = convert(ref, start, self.tsData, cs, seq)

        flag = 0
        if strand_g == False:
            flag = 16

        a = pysam.AlignedSegment()
        a.query_name = read.query_name
        a.flag = flag
        a.reference_id = self.reflist.index(chrom_g)
        a.reference_start = r_st_g
        a.mapping_quality = 20
        a.cigarstring = cigar_g
        a.query_sequence = seq_g
        # set referencename
        a.set_tag("XS", ref)

        for tag, value in read.get_tags():

            if containHardClip:
                if tag == "MM" or tag == "MM":
                    print("MM tag removed ts")
                    a.set_tag("XR", 1)
                    continue

            a.set_tag(tag, value)

        return a


def remove(read, tag_to_remove):
    if read.has_tag(tag_to_remove):
        read.set_tag(tag_to_remove, None, value_type=None)


def checkOverlap(read, chrom, strand, r_st, r_en):
    if read.reference_name != chrom:
        return False
    if (not read.is_reverse) != strand:
        return False
    #
    start = read.reference_start
    end = read.reference_end
    return (start <= r_en) and (r_st <= end)


def genomicProcess(read, aligner, in_bamfile):
    resultExsist = False
    for hit in aligner.map(read.seq):

        # take tophit only
        chrom = hit.ctg
        strand = hit.strand
        flg = 0
        if not strand:
            flg = 16

        r_st = hit.r_st
        r_en = hit.r_en
        q_st = hit.q_st
        q_en = hit.q_en
        #
        readlen = len(read.seq)
        #
        if strand:
            pre = q_st
            post = readlen - q_en
        else:
            pre = readlen - q_en
            post = q_st

        cigarStr = hit.cigar_str
        if pre > 0:
            cigarStr = str(pre) + "S" + cigarStr
        if post > 0:
            cigarStr = cigarStr + str(post) + "S"

        read.flag = flg
        read.reference_id = in_bamfile.get_tid(hit.ctg)
        read.reference_start = r_st
        read.mapping_quality = hit.mapq
        read.cigarstring = cigarStr
        read.set_tag('NM', hit.NM)

        overlap = checkOverlap(read, chrom, strand, r_st, r_en)
        if not overlap:
            remove(read, "MM")
            remove(read, "ML")
            print("MM tag removed map")
            read.set_tag("XR", 1)

        resultExsist = True
        break
    #
    if (resultExsist == False):

        cigar_tuples = read.cigar
        modified_cigar = []
        containHardClip = False
        for op, length in cigar_tuples:

            if op == 5:  # hard clip
                containHardClip = True
            if op == 2 and length >= 80:  # 2: deletion
                modified_cigar.append((3, length))  # 3: N (skip)
            else:
                modified_cigar.append((op, length))

        if containHardClip:
            remove(read, "MM")
            remove(read, "ML")
            print("MM tag removed ge")
            read.set_tag("XR", 2)

        read.cigar = modified_cigar

    return read


import os
import time
import pysam


def is_bai_outdated(bam_file, bai_file):
    """Check if the BAI file is older than the BAM file or if it doesn't exist."""
    if not os.path.exists(bam_file) or not os.path.exists(bai_file):
        return True  # The BAI file needs updating if it doesn't exist.
    bam_mtime = os.path.getmtime(bam_file)
    bai_mtime = os.path.getmtime(bai_file)
    return bai_mtime < bam_mtime


def check_and_update_bam_index(bam_file, sleeptime=60, max_tries=30):
    bai_file = bam_file + ".bai"
    """Check and potentially update the BAM file's index, with a maximum number of attempts."""
    for _ in range(max_tries):
        if is_bai_outdated(bam_file, bai_file):
            # Update the BAM file index using pysam.
            pysam.index(bam_file)

        # Wait for 1 minute.
        time.sleep(sleeptime)
        # Exit the loop early if the BAI file is no longer outdated.
        if not is_bai_outdated(bam_file, bai_file):
            break


import subprocess
import os
import statistics
import numpy as np


def run_convert(inbam, outbam, ref, convertMatrix):
    aligner = mp.Aligner(ref, preset="splice", k=14, extra_flags=0x100, best_n=1)
    sqlist = []
    reflist = []
    for record in SeqIO.parse(ref, 'fasta'):
        sqlist.append({'LN': len(record), 'SN': record.id})
        reflist.append(record.id)
    header = {'HD': {'VN': '1.0'},
              'SQ': sqlist}

    mapped_reads_count = 0
    unmapped_reads_count = 0
    read_lengths = []

    in_bamfile = pysam.AlignmentFile(inbam, "rb")
    out_bamfile = pysam.AlignmentFile(outbam, "wb", header=header)
    converter = Converter(reflist, convertMatrix)
    chromosome_counts = {}

    n = 0
    m = 0
    gm = 0
    for read in in_bamfile:

        # for counter
        if read.is_unmapped:
            unmapped_reads_count += 1
        else:
            mapped_reads_count += 1
            read_lengths.append(read.query_length)

        if not read.is_unmapped and read.seq is not None:

            mappedref = read.reference_name
            if mappedref in reflist:
                #    read mapped to genome no conversion
                c_read = genomicProcess(read, aligner, in_bamfile)
                gm += 1
            else:
                #    read mapped to transcript convert to genomic coordinate
                c_read = converter.convert_read(read)
                m += 1

            out_bamfile.write(c_read)

            # for reads counter
            chromosome = None
            rid = c_read.reference_id
            if (rid >= 0) and (rid < len(reflist)):
                chromosome = reflist[rid]
            if chromosome is not None:
                if chromosome in chromosome_counts:
                    chromosome_counts[chromosome] += 1
                else:
                    chromosome_counts[chromosome] = 1

        n += 1

        if n % 1000 == 0:
            print(m, gm, n)

    in_bamfile.close()
    out_bamfile.close()

    sorted_bam = os.path.splitext(outbam)[0] + "_sorted.bam"
    pysam.sort("-o", sorted_bam, outbam)
    check_and_update_bam_index(sorted_bam)

    # for read counter
    if len(read_lengths) > 0:
        median_length = statistics.median(read_lengths)
        min_length = min(read_lengths)
        max_length = max(read_lengths)
    else:
        median_length = min_length = max_length = 0

    read_lengths_np = np.array(read_lengths)
    maxv = max(max(read_lengths), 50000)
    hist, bin_edges = np.histogram(read_lengths_np, bins=range(0, maxv + 101, 100))
    with open(outbam[0:-3] + "_reads_summary.txt", "w") as f:

        f.write("#Reads stats\n")
        f.write("Mapped Reads: {}\n".format(mapped_reads_count))
        f.write("Unmapped Reads: {}\n".format(unmapped_reads_count))
        f.write("Median Length: {}\n".format(median_length))
        f.write("Minimum Length: {}\n".format(min_length))
        f.write("Maximum Length: {}\n".format(max_length))

        f.write("\n")
        f.write("#Bin_Start Bin_End Count\n")
        for i in range(len(hist)):
            bin_start = bin_edges[i]
            bin_end = bin_edges[i + 1] - 1
            count = hist[i]
            f.write(f"{bin_start} {bin_end} {count}\n")

        f.write("\n")
        f.write("#Chromosome counts\n")
        chrlist = sorted(list(chromosome_counts.keys()))
        for chromosome in chrlist:
            count = chromosome_counts[chromosome]
            f.write(f"{chromosome}\t{count}\n")


import os


def test():
    ref = "/share/reference/hg38.fa"

    # inbam = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/U87transcript.bam"
    # outbam = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/U87_ge.bam"
    # convertMatrix = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/bam_pass/U87wt.matrix"
    #
    # run_convert(inbam,outbam,ref,convertMatrix)
    #
    # inbam = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/U87ivttranscript.bam"
    # outbam = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/U87ivt_ge.bam"
    # convertMatrix = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/bam_pass/U87ivt.matrix"
    #
    # run_convert(inbam,outbam,ref,convertMatrix)

    inbam = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/U87_STM2457ivttranscript.bam"
    outbam = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/U87_STM2457ge.bam"
    convertMatrix = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/bam_pass/U87_STM2457.matrix"

    run_convert(inbam, outbam, ref, convertMatrix)

# test_stats()

