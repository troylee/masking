#!/usr/bin/env python
import string, sys
from scipy.io import savemat

scplst=[]
fin=open('/home/li-bo/aurora2/data/lists/wav2mat_train_all.scp')
while True:
	sr=fin.readline()
	if sr=='':break
	sr=sr.strip()
	lst=sr.split(' ')
	sr=lst[2]
	fname=sr[sr.rfind('/')+1:] # remove path
	fname=fname[:fname.rfind('.')] # remove suffix
	if fname.startswith('clean_'):
		scplst.append(fname)
fin.close()

assert(len(scplst)==8440)

alignlst=[{} for i in range(8440)]
fin=open('/home/li-bo/aurora2/exps/cleantr/a00_hmm/mono1b/pdf_align/train.pdf')
cnt=0
while True:
	sr=fin.readline()
	if sr=='':break
	sr=sr.strip()
	lst=sr.split(' ')
	fname=lst[0]
	statelst=[]
	for i in range(1, len(lst)):
		statelst.append(int(lst[i]))
	idx=scplst.index(fname)
	alignlst[idx]['name']=fname
	alignlst[idx]['states']=statelst
	cnt+=1
fin.close()
assert(cnt==8440)

res={}
res['pdfaligns']=alignlst
savemat('/home/li-bo/aurora2/exps/cleantr/a00_hmm/mono1b/pdf_align/train_pdf.mat', res)

