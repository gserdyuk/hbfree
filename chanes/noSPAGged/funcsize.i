
C values fo internal array sizes:
C MAXU_MD  - max number of voltages per model
C MAXD_MD  - max number of voltage derivetives per model
C MAXKN - max number of frequencies ( almost all 20-th in this program
C         is MAXKN
C MAXKN1 - max size of extended frequency grid (for derivatives)
C size W in common/blw1/ is MAXKN, size W1 - MAXKN1
       integer MAXU_MD, MAXD_MD, MAXKN, MAXKN1
       parameter (MAXU_MD = 10)
       parameter (MAXD_MD = 10)
       parameter (MAXKN = 20)
       parameter (MAXKN1 = 200)


C B1_SIZE - descriptins of B1 and B2 are in FTMAS2, size of B1 = max size of FFT
C B2_SIZE - size of B2 - shoudl be enough to store sparce 2-dim vector
C           (matrix). max - B1_SIZE*B1_SIZE (dence case)

       integer B1_SIZE, B2_SIZE
       parameter (B1_SIZE = 256)
       parameter (B2_SIZE = 256*16)

