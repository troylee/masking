#!/usr/bin/env python
import string, sys

if len(sys.argv)!=3:
	print 'Usage: script <wavscp> <resscp>'
	sys.exit(1)
	
fin=open(sys.argv[1])
fout=open(sys.argv[2],'w')
while True:
	sr=fin.readline()
	if sr=='':break
	sr=sr.strip()
	res=string.replace(sr, '/home/li-bo/aurora2/data/wav/', '');
	
	# clean wav filename
	clr=''
	lst=res.split('/')
	if 'train' in sr:
		clr='/home/li-bo/aurora2/data/wav/'+lst[0]+'/clean/'+lst[2]
	else:
		if lst[1].startswith('clean'):
			clr=sr
		else:
			nid=lst[1][1]
			clr='/home/li-bo/aurora2/data/wav/'+lst[0]+'/clean'+nid+'/'+lst[2]

	# result txt filename
	res=string.replace(res, '.wav', '.txt')
	sid=res.rfind('/')
	print >>fout, sr, clr, res[:sid]+'_'+res[sid+1:]
fin.close()
fout.close()

