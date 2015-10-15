
import numpy as np
cimport numpy as np
import struct
import sys

# maps plink binary represntation of genotypes to an unsigned integer
# missing values are coded by `3`.
cdef dict genomap = dict([('00',0),('01',1),('11',2),('10',3)])

def load(file):

    cdef int n, l, i, Nindiv, Nsnp, Nbytes
    cdef bytes line
    cdef str checkA, checkB, checkC, bytestr
    cdef np.ndarray genotype

    # number of individuals
    handle = open(file+'.fam','r')
    for i, li in enumerate(handle):
        pass
    Nindiv = i+1

    # Number of bytes to read in at a time
    Nbytes = Nindiv/4+(Nindiv%4>0)*1

    # number of SNPs
    handle = open(file+'.bim','r')
    for i, li in enumerate(handle):
        pass
    Nsnp = i+1

    tobit = lambda x: ''.join([bin(i)[2:].zfill(8)[::-1] for i in struct.unpack('<%sB'%Nbytes, x)])
    genotype = np.zeros((Nindiv,Nsnp),dtype='uint8')

    # open the bed file
    handle = open(file+'.bed','rb')

    # check if the file is a valid plink bed file
    line = handle.read(1)
    checkA = bin(struct.unpack('<B', line)[0])[2:].zfill(8)[::-1]
    line = handle.read(1)
    checkB = bin(struct.unpack('<B', line)[0])[2:].zfill(8)[::-1]
    line = handle.read(1)
    checkC = bin(struct.unpack('<B', line)[0])[2:].zfill(8)[::-1]

    if checkA!="00110110" or checkB!="11011000":
        print "This is not a valid bed file"
        handle.close()
        sys.exit(2)

    # parse the bed file        
    for l from 0 <= l < Nsnp:
        line = handle.read(Nbytes)
        bytestr = tobit(line)
        for n from 0 <= n < Nindiv:
            genotype[n,l] = genomap[bytestr[2*n:2*n+2]]

    handle.close()

    return genotype

