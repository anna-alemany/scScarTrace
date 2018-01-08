import sys, os
from pandas.io.parsers import read_csv

dfname = []
df = []
dflabel = []
try:
    num = int(sys.argv[1])
    if num <=1:
        sys.exit('What do you want to merge?')
    for i in range(num):
        dfname.append(sys.argv[2*i+2])
        print i, dfname[-1]
        df.append(read_csv(sys.argv[2*i+2], sep = '\t', index_col = 0))
        dflabel.append(sys.argv[2*i+3])
    how2merge = sys.argv[-2]
    out = sys.argv[-1]
except:
    print sys.argv[1:]
    sys.exit('Please, (1) give number of input files; (2) input df to merge; (3) how to merge; (4) output root')

print 'Read all'

for i in range(num):
    df[i].columns = [col + '-' + dflabel[i] for col in df[i].columns]

pm = df[0].merge(df[1], how = how2merge, left_index = True, right_index = True)
pm = pm.fillna(0)
if num > 2:
    for i in range(2, num):
        pm = pm.merge(df[i], how = how2merge, left_index = True, right_index = True)
        pm = pm.fillna(0)

pm = pm.loc[pm.sum(axis=1).sort_values(ascending=False).index]

pm.to_csv(out + '.tsv', sep = '\t')
f = open(out + '.log', 'w')
print >> f, dfname
print >> f, dflabel
f.close()

 
