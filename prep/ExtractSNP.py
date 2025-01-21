import gzip

def cond(ref,alt):

    shift = 0
    cond1 = (ref == "T" and alt =="C")
    cond2 = (ref == "A" and alt == "G")
    if cond1 or cond2:
        return True,shift

    if "," in alt:
        alt = ref.split(",")
        cond1 = (ref == "T" and "C" in alt)
        cond2 = (ref == "A" and "G" in alt)
    if cond1 or cond2:
        return True,shift

    cond1 = (ref == "T" and alt =="TC")
    cond2 = False
    if len(alt) ==2:
        if alt[1] == "G":
            cond2 = True,shift
            shift = 1
    if cond1 or cond2:
        return True,shift

    cond1 = (ref == "C" and alt == "T")
    cond2 = (ref == "G" and alt == "A")
    if cond1 or cond2:
        return True, shift

    return False,shift


import pandas as pd
def prepareDB(invcf,outdir):

    tclist = []
    aglist = []
    with gzip.open(invcf, 'rt') as f:
        n=0
        lastchr = None
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            data =  line.split("\t")
            chr = "chr"+data[0]
            pos = int(data[1])
            id = data[2]
            ref = data[3]
            alt = data[4]
            condition,shift = cond(ref, alt)
            pos = pos + shift
            if condition:
                print(chr, pos, id, ref, alt)
                if ((ref == "T") and (shift == 0)) or (ref == "C") :
                    tclist.append((pos,id))
                else:
                    aglist.append((pos, id))
                n+=1

                if lastchr is not None and (lastchr!=chr):
                    df = pd.DataFrame(tclist, columns=['pos', 'snpId'])
                    df.to_parquet(outdir + '/'+lastchr+'_tc.parquet')
                    df = pd.DataFrame(aglist, columns=['pos', 'snpId'])
                    df.to_parquet(outdir + '/'+lastchr+'_ag.parquet')
                    tclist = []
                    aglist = []
                lastchr = chr

            # if n== 1000:
            #     break

        df = pd.DataFrame(tclist, columns=['key', 'snpId'])
        df.to_parquet(outdir + '/' + lastchr + '_tc.parquet')
        df = pd.DataFrame(aglist, columns=['key', 'snpId'])
        df.to_parquet(outdir + '/' + lastchr + '_ag.parquet')

def test():

    invcf = '/mnt/share/ueda/RNA004/resource/00-All.vcf.gz'
    outdir = "/mnt/share/ueda/RNA004/resource/dbSNP151TConly/"
    prepareDB(invcf,outdir)

#test()cd /cp ./