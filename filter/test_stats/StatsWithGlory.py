import mappy as mp
import pandas as pd

# Function to read data from a file
def stats(vcf1,knownPos):

    bed_df = pd.read_csv(vcf1, header=None, sep='\t')
    counter = {}
    counterknown = {}
    cnt = 0
    drachcnt = 0
    for index, row in bed_df.iterrows():

        chr = row[0]
        pos = row[1]
        alt = row[4]
        filter = row[6]
        info = row[7]
        if "rescued=True" in info and "knownSites=True" in info:
            continue

        key0 = str(chr) + ":" + str(pos)
        inKnownDB = key0 in knownPos
        drach = False
        if alt != "a":
            continue
        if alt == "a" and "knownMotif=True" in info:
            pass
        else:
            continue

        if filter:

            if inKnownDB:
                if alt in counterknown:
                    counterknown[alt] = counterknown[alt] + 1
                else:
                    counterknown[alt] = 1

            if alt in counter:
                counter[alt] = counter[alt] + 1
            else:
                counter[alt] = 1
            cnt += 1




        if cnt % 1000 == 0:
            print(cnt, drachcnt,counter,counterknown)
    print(drachcnt)
    print(counter)
    print(counterknown)

# File paths
file1_path = '/mnt/share/ueda/RNA004/resource/glory.bed'
# file2_path = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"
file2_path = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v02/HEK293T_DR13/HEK293T_DR13/filter_result.vcf"
# file4_path = '/mnt/share/ueda/RNA004/nanoZero/unfilter_result.vcf'
knownPos = {}
bed_df = pd.read_csv(file1_path, header=None, sep='\t')
cnt = 0
print(bed_df)
for index, row in bed_df.iterrows():

    chr = row[0]
    pos = row[2]
    key = str(chr) + ":" + str(pos)
    print(key)
    knownPos[key] = 1
    cnt+=1

stats(file2_path,knownPos)
print(cnt)





