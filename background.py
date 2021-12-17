from collections import Counter
import numpy as np
conv = {"A":0, "C":1, "G":2, "T":3}

#Cleans and truncates Sequeces to be of equal length
def clean_strings(Seqs):
    rmv =  [i.replace("... <Preview truncated at 128 characters>","") for i in Seqs ]
    minSq = min([len(i) for i in Seqs])
    return [i[:minSq] for i in rmv]

#Calculates nucletide frequencies for A,C,G,T    
def zero_order_freq(ls):
    x = "".join(ls)
    tmp = dict(Counter(x))
    nucls = len(x)
    return np.array([tmp["A"]/nucls,tmp["C"]/nucls,tmp["G"]/nucls,tmp["T"]/nucls])

#Calculates count of 2-mers in a sequence : AA,AC,AG,AT,CA .... 
def get_first_order_count(Seqs):
    dct = {}
    def add_seq(ls):
        for i in range(len(ls)-1):
            kmr = ls[i:i+2]
            if dct.__contains__(kmr):
                dct[kmr] += 1
            else:
                dct[kmr] = 1
    for seq in Seqs:
        add_seq(seq)
    return dct
#Return  2-mer frequency matrix   
def first_order_freq(Seqs):
    count_dc = get_first_order_count(Seqs)
    freq = np.zeros((4,4))
    ttl = sum(count_dc.values())
    for key in count_dc:
        freq[conv[key[0]]][conv[key[1]]] = count_dc[key] / ttl
    return freq
#Given a frequency matrix of some order score the probability of finding
# each nucleotide at each position in a sequence    
def create_nucl_background(Seqs,order,freq_m):
    M,N  = len(Seqs),len(Seqs[0])
    if order == 0:
        bgm = np.zeros((M,N))
        for i in range(M):
            for j in range(N):
                bgm[i][j] = freq_m[conv[Seqs[i][j]]]
        return bgm
    elif order == 1:
        bgm = np.zeros((M,N))
        avg_vals = {"A":freq_m[:,conv["A"]].mean(),"C":freq_m[:,conv["C"]].mean(),"G":freq_m[:,conv["G"]].mean(),"T":freq_m[:,conv["T"]].mean()}
        for i in range(M):
            for j in range(N):
                if j == 0:
                    bgm[i][j] = avg_vals[Seqs[i][j]]
                else:
                    bgm[i][j] = freq_m[conv[Seqs[i][j-1]]][conv[Seqs[i][j]]]
        return bgm
    else:
        return ValueError("Order must be 0 or 1. Other orders not supported")
#Calculates probability of observing each kmer given the background per nucleotide        
def create_kmer_background(bgm,kmr:int):
    bkm = np.zeros((bgm.shape[0],bgm.shape[1]-kmr))
    for i in range(bkm.shape[0]):
        for j in range(bkm.shape[1]):
            bkm[i][j] = bgm[i][j:j+kmr].prod()
    return bkm
def get_background(Seqs,kmr,freq_m,order):
    bgm = create_nucl_background(Seqs,order,freq_m)
    bkm = create_kmer_background(bgm,kmr)
    return bkm
