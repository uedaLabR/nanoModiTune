def insertN(cigar,tsName,start,tsData,debug=False):

    # findout N position as list
    nlist = []
    tsInfo = tsData[tsName]
    exons = tsInfo.exons
    exonLen = len(exons)
    nlensum = 0
    for n in range(1,exonLen):

       s0,e0 = exons[n-1]
       s1, e1 = exons[n]
       pos = e0-start-nlensum
       nlen = s1-e0
       if pos > 0:
           nlist.append((pos,nlen))
           nlensum = nlensum + nlen


    return _insertN(cigar,nlist,debug)

def _insertN(cigar,nlist,debug):

    # Insert N
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    cigarlist = []
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator != 3:
            cigarlist.append((cigaroprator, cigarlen))

    for pos,nlen in nlist:
        if debug:
            print("insert N in",pos,nlen)
            print(toCigarStr(cigarlist))
        cigarlist = insertN_each(pos,nlen,cigarlist)

    return toCigarStr(cigarlist)

def toCigarStr(cigarlist):

    cgstr =""
    for cigaroprator, cigarlen in cigarlist:

        if cigaroprator == 3:  # N
            cgstr = cgstr + str(cigarlen) +"N"
        elif cigaroprator == 0 :  # match M
            cgstr = cgstr + str(cigarlen) + "M"
        elif cigaroprator == 2:  # Del
            cgstr = cgstr + str(cigarlen) + "D"
        elif cigaroprator == 1:  # Ins
            cgstr = cgstr + str(cigarlen) + "I"
    return cgstr

def insertN_each(pos,nlen,cigarlist):

    newlist = []
    refpos = 0
    cigaropratorN = 3
    beforeadd = True
    for cigaroprator, cigarlen in cigarlist:

        if cigaroprator == 3:  # N

            newlist.append((cigaroprator, cigarlen))

        elif cigaroprator == 0 or cigaroprator == 2:   # match or Del

            if  (refpos <= pos) and (pos <= refpos + cigarlen) and beforeadd:

                left = pos-refpos
                if left > 0:
                    newlist.append((cigaroprator, left))
                newlist.append((cigaropratorN, nlen))
                #print("insertn",pos, nlen)

                right = cigarlen-left
                if right > 0:
                    newlist.append((cigaroprator, right))
                beforeadd = False
            else:
                newlist.append((cigaroprator, cigarlen))

            refpos = refpos + cigarlen



        elif cigaroprator == 1:  # Ins

            newlist.append((cigaroprator, cigarlen))

    return newlist
