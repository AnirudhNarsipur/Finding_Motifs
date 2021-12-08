import random 
import numpy as np
from numpy.core.defchararray import count
conv = {"A":0, "C":1, "G":2, "T":3}
def score_motif(pr_m,seq):
    s = 1
    for i in range(0,len(seq)):
        s*=pr_m[conv[seq[i]]][i]
    return s
def GibbsSampler(Seqs : "list[str]", k:int, numSeq:int, seqLen:int,itr:int):
    rand_starts = random.choices(range(0,seqLen-k+1),k=numSeq)
    choose_random_seq = lambda : random.randint(0, numSeq - 1)
    count_matrix = np.zeros((4,k))
    profile_matrix = np.zeros((4,k))
    last_k_pos = seqLen - k 
    def generate_count_matrix(exclude):
        for i in range(0,numSeq):
            if i == exclude:
                continue
            else:
                kstart = rand_starts[i]
                for j in range(kstart,kstart+k):
                    count_matrix[conv[Seqs[i][j]]][j-kstart] += 1
    def set_excluded_motif(exclseq : int):
        present_kmers = np.zeros((last_k_pos + 1))
        for i in range(0,last_k_pos):
            present_kmers[i] = score_motif(profile_matrix,Seqs[exclseq][i:i+k])
        sm = sum(present_kmers)
        present_kmers=present_kmers/sm
        rand_starts[exclseq] = random.choices(range(0,last_k_pos+1),present_kmers,k=1)[0]
    exclude_seqs = random.choices(range(0,numSeq),k=itr)
    for i in range(itr):
        exclude_seq : int = exclude_seqs[i]
        generate_count_matrix(exclude_seq)
        profile_matrix = (count_matrix + 1)/(4+(numSeq-1))
        set_excluded_motif(exclude_seq)
        count_matrix = np.zeros((4,k))
    return profile_matrix
if __name__ == "__main__": 
    seqs = ["ATCGAA","ACTTCG"]
    print(GibbsSampler(seqs,k=3,numSeq=len(seqs),seqLen=6,itr=3000))
#Rename column in Pandas dataframe df
#s