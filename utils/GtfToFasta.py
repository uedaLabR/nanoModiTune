import mappy as mp
from Bio.SeqFeature import SeqFeature, FeatureLocation

def toSt(locations):

    s = ""
    idx = 0
    for location in locations:

        if idx > 0:
            s = s+","
        s = s+ "{0}:{1}".format(location.start-1, location.end-1)

        idx+=1

    return s

def formatData(d):

    transcript_id, chr,strand, locations= d
    poss = toSt(locations)
    line = "{0}\t{1}\t{2}\t{3}".format(transcript_id,chr, strand, poss)
    return line

def toseq(a,chr,strand,locations):

    ss = ""
    for location in locations:

        s = a.seq(chr, location.start-1, location.end)
        if s is not None:
            ss = ss + s

    if strand == -1:
        ss = mp.revcomp(ss)

    return ss

def checkConsist(locations):

    first = locations[0]
    last = locations[len(locations)-1]
    # print("b4",locations)
    if first.start > last.start:

        locations.reverse()

    if len(locations) > 1:
        first = locations[0]
        second = locations[1]
        # if inconsistant remove first exon
        if first.end > second.start:
            print("something went wrong")
            locations.pop(0)

    return locations

# def extends3UTR(locations,strand):
#
#     first = locations[0]
#     last = locations[len(locations)-1]
#     margin = 500
#     if strand:
#         location = FeatureLocation(last.start, last.end + margin, strand=strand)
#         locations.pop()
#         locations.insert(location)
#
#     else:
#
#         location = FeatureLocation(first.start-margin, first.end, strand = strand)
#         locations.pop(0)
#         locations.insert(location)
#
#     return locations

def getID(data):

    info = data[8].split(";")
    if "reference_id" in data[8]:
        for i in info:
            if "reference_id" in i:
                return i.replace("reference_id", "").replace('"', '').replace(' ', '')

    for i in info:
        if "ID=" in i:
            return i.replace("ID=","")
        if "transcript_id" in i:
            return i.replace("transcript_id", "").replace('"','').replace(' ','')
    return None

def getTranscript(gtf,excludeKnown):

    rdata = []
    with open(gtf) as handle:

        idbefore = None
        locations = []
        transcript_id = None

        transcript_idset = set()

        for s_line in handle:

            if s_line.startswith("#"):
                continue


            if excludeKnown:
                if "reference_id" in s_line:
                    continue

            data = s_line.split("\t")
            if len(data) < 3:
                break

            type = data[2]
            if type == "gene":

                chr = data[0]
                strand = 1
                if data[6] == "-":
                    strand = -1
                transcript_id = getID(data)
                if (transcript_id is not None) and (transcript_id not in transcript_idset):
                    if len(locations) > 0:
                        checkConsist(locations)
                        transcript_id_b4, chr_b4, strand_b4 = idbefore
                        rdata.append((transcript_id_b4, chr_b4, strand_b4, locations))
                        transcript_idset.add(transcript_id)

                idbefore = transcript_id, chr, strand
                locations = []

            elif type == "transcript":

                chr = data[0]
                strand = 1
                if data[6] == "-":
                    strand = -1
                transcript_id = getID(data)
                if (transcript_id is not None) and (transcript_id not in transcript_idset):
                    if len(locations) > 0:
                        checkConsist(locations)
                        transcript_id_b4, chr_b4, strand_b4 = idbefore
                        rdata.append((transcript_id_b4, chr_b4, strand_b4, locations))
                        transcript_idset.add(transcript_id)

                idbefore = transcript_id, chr, strand
                locations = []


            elif type == "exon":

                location = FeatureLocation(int(data[3]), int(data[4]), strand=strand)
                locations.append(location)
                strand = location.strand

        if (transcript_id is not None) and (transcript_id not in transcript_idset):
            if len(locations) > 0:
                checkConsist(locations)
                transcript_id_b4, chr_b4, strand_b4 = idbefore
                rdata.append((transcript_id_b4, chr_b4, strand_b4, locations))
                transcript_idset.add(transcript_id)

    return rdata

def copy_large_file(src_file_path, dst_file_path, buffer_size=1024*1024):
    """
    Function to copy a large file in chunks.

    :param src_file_path: Path of the source file to be copied
    :param dst_file_path: Path of the destination file
    :param buffer_size: Size of the data to read at once (in bytes)
    """
    with open(src_file_path, 'rb') as src_file:
        with open(dst_file_path, 'wb') as dst_file:
            while True:
                chunk = src_file.read(buffer_size)
                if not chunk:
                    break
                dst_file.write(chunk)


def register(checkdict,ts):

    transcript_id, chr, strand, locations = ts
    chrstrandkey = chr+":"+str(strand)
    if chrstrandkey in checkdict:
        tlist_d = checkdict[chrstrandkey]
    else:
        tlist_d = {}
        checkdict[chrstrandkey] = tlist_d

    first = locations[0]
    last = locations[len(locations)-1]
    s = first.start
    e = last.end
    binkey = (s+e)//10000

    if binkey in tlist_d:
        tlist = tlist_d[binkey]
    else:
        tlist = []
        tlist_d[binkey] = tlist

    tlist.append(((s,e),locations))

def equalwithmajin(lc1,lc2,margin):

    b1 =  abs(lc1.start-lc2.start)<=margin
    b2 = abs(lc1.end - lc2.end) <= margin
    return b1 and b2

def checkExact(locations1,locations2):

    len1 = len(locations1)
    len2 = len(locations2)
    if len1 != len2:
        return False
    else:

        if len1 == 1:

            margin = 30
            return equalwithmajin(locations1[0], locations2[0], margin)

        else:

            margin = 5
            for n in range(len1):

                lc1 = locations1[n]
                lc2 = locations2[n]
                if not equalwithmajin(lc1,lc2,margin):
                    return False

    return True

def isExsist(tlist_d,tp):

    (s, e), locations = tp
    binkey = (s+e)//10000
    if binkey in tlist_d:

        tlist =  tlist_d[binkey]
        for tp1 in tlist:
            (s1, e1), locations1 = tp1
            #overlap
            if s<=e1 and s1 <=e:
                exactmatch = checkExact(locations,locations1)
                if exactmatch:
                    return True
    else:
        return False


def nonoverrap(checkdict,ts):

    transcript_id, chr, strand, locations = ts
    first = locations[0]
    last = locations[len(locations)-1]
    s = first.start
    e = last.end
    chrstrandkey = chr+":"+str(strand)
    if chrstrandkey in checkdict:

        tlist_d = checkdict[chrstrandkey]
        exsist = isExsist(tlist_d,((s,e),locations))
        if exsist:
            return False

    return True


def mergeTranscript(stingtiedata):

    rdata = []
    firstlist = stingtiedata[0]
    rdata.extend(firstlist)
    checkdict = {}
    for ts in firstlist:
        register(checkdict,ts)

    #
    cntoverlap = 0
    for n in range(1,len(stingtiedata)):

        data = stingtiedata[n]
        for ts in data:

            if nonoverrap(checkdict,ts):

                transcript_id, chr, strand, locations = ts
                ts = transcript_id+"_"+str(n), chr, strand, locations
                rdata.append(ts)
                register(checkdict, ts)
            else:
                cntoverlap+=1
        print(len(data),cntoverlap)
    return rdata


import os
import glob
def process(gtf,gtf_stringtie,ref,outmatrix,outfasta):

    print("start process")
    print("ref", ref)

    print("reading gtf1",gtf)
    data = getTranscript(gtf,False)

    is_directory = os.path.isdir(gtf_stringtie)
    if is_directory:

        files = glob.glob(gtf_stringtie+"/*.gtf")
        stingtiedata = []
        if len(files) > 1:
            for file in files:
                print("reading gtf2", file)
                data2 = getTranscript(file, True)
                stingtiedata.append(data2)
            data_add = mergeTranscript(stingtiedata)
        else:
            data_add = getTranscript(files[0], True)
    else:
        data_add = getTranscript(gtf_stringtie, True)

    data.extend(data_add)

    print("end reading file")
    print("write matrix")
    with open(outmatrix, mode='w') as f:

        for d in data:
            s = formatData(d)
            f.write(s + "\n")

    if os.path.exists(outfasta):
        os.remove(outfasta)


    print("write sequence")
    a = mp.Aligner(ref)
    with open(outfasta, mode='a') as f:

        for transcript_id, chr, strand, locations in data:
            seq = toseq(a, chr, strand, locations)
            f.write(">" + transcript_id + "\n")
            f.write(seq + "\n")


def run():

    gtf =  "/mnt/share/ueda/RNA004/pipeline/gencode.v44.annotation.gff3"
    ref = "/share/reference/hg38.fa"

    gtf_stringtiedir = "/mnt/share/ueda/RNA004/U87/U87_result/refs"
    outmatrix = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/bam_pass/U87wt.matrix"
    outfasta = "/mnt/share/ueda/RNA004/U87/U87_WT/20231215_1505_MN32625_FAX70236_31c89747/bam_pass/U87wt.fa"
    process(gtf,gtf_stringtiedir, ref, outmatrix, outfasta)

    # gtf_stringtie = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/bam_pass/U87ivt.gtf"
    # outmatrix = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/bam_pass/U87ivt.matrix"
    # outfasta = "/mnt/share/ueda/RNA004/U87/U87_IVT/20231227_1535_MN32625_FAX73794_2cf3868f/bam_pass/U87ivt.fa"
    # process(gtf,gtf_stringtie,ref,outmatrix,outfasta)
    #
    # gtf_stringtie = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/bam_pass/U87_STM2457.gtf"
    # outmatrix = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/bam_pass/U87_STM2457.matrix"
    # outfasta = "/mnt/share/ueda/RNA004/U87/U87_STM2457/U87_STM2457/20240216_1602_MN32625_FAX70218_f0b7f679/bam_pass/U87_STM2457.fa"
    # process(gtf,gtf_stringtie,ref,outmatrix,outfasta)

# run()