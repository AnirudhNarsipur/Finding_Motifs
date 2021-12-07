import numpy as np
comp = {"A":["U"],"U":["A","G"],"C":["G"],"G":["C","U"]}
def is_comp(nucl,nucl2):
    return nucl2 in comp[nucl]

def get_rev(seq,i):
    return seq[-i-1]
def al(seq):
    sqlen = len(seq)
    mtrx = np.zeros((sqlen,sqlen)).astype(int)
    for i in range(sqlen):
        rev = seq[-i-1]
        if is_comp(seq[0] ,rev):
            mtrx[i][0] = 1
        if is_comp(seq[-1],seq[i]):
            mtrx[0][i] = 1
    for i in range(1,sqlen):
        for j in range(1,sqlen):
            if is_comp(seq[-i-1], seq[j]):
                mtrx[i][j] = mtrx[i-1][j-1] + 1 
    return mtrx
    
def print_mtrx(seq,mtrx):
    print(" ",end=" ")
    for i in seq:
        print(i,end=" ")
    print()
    rev = seq[::-1]
    for i in range(len(mtrx)):
        print(rev[i],end=" ")
        for j in range(len(mtrx[0])):
            print(mtrx[i][j],end=" ")
        print()

if __name__ == "__main__":
    sq = "CAGAUUUACUAGUACGUAAUUG"
    print_mtrx(sq,al(sq))

    