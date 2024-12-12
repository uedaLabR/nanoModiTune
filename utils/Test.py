file_path = "/mnt/ssdnas07/nanozero/rna/nanomoditune_v01/HEK293T_DR13/HEK293T_DR13/unfilter_result.vcf"
file_path = "/mnt/share/ueda/RNA004/hek293/result_filter.vcf"
cnt=0
with open(file_path, "r") as file:
    for line in file:
        if "17596" in line:
            print(line.strip())
            cnt+=1
print(cnt)
