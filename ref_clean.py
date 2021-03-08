import sys
import csv
import pandas as pd

inf = pd.read_csv(sys.argv[1],sep='\t')
dmr = pd.read_csv(sys.argv[2],sep='\t')
dict={}
ref_chr=pd.to_numeric(inf['chr'].iloc[0])
ref_start=pd.to_numeric(inf['start'].iloc[0])
ref_end=pd.to_numeric(inf['end'].iloc[0])
chrstart=str(ref_chr)+'_'+str(ref_start)
dict[chrstart]=ref_end

for index,row in inf.iterrows():
	chr=pd.to_numeric(row['Chr'])
	start=pd.to_numeric(row['Start'])
	end=pd.to_numeric(row['End'])
	if ref_chr==chr:
		if start <= ref_end:
			continue
		elif start > ref_end:
			ref_start=start
			ref_end=end
			chrstart=str(ref_chr)+'_'+str(ref_start)
			dict[chrstart]=ref_end
	else:
		ref_chr=chr
		ref_start=start
		ref_end=end
		chrstart=str(ref_chr)+'_'+str(ref_start)
		dict[chrstart]=ref_end
for key in dict.keys():
	print str(key)+"\t"+str(dict[key])
#df = pd.DataFrame(list(dict.items()),columns = ['Chrstart','End'])
#df[['Chr','Start']]=df.Chrstart.str.split(expand=True)
