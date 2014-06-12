path=/home/li-bo/aurora2/data/noisemask/ftrs
curpath=`pwd`

for ftr in IBMFBankSpeech TBMFBankSpeech
do
	# training
	for data in cleantr #multitr
	do
		cat ../fbank_e_d_a_kaldi/${data}.scp | awk '{print $1}' > tmp.scp
		copy-feats-from-text --dim=72 --data-suffix=txt --data-directory=ftrs/${ftr}/train tmp.scp ark,scp:${path}/${ftr}_kaldi/${data}.ark,${path}/${ftr}_kaldi/${data}.scp
		
		compute-cmvn-stats scp:${path}/${ftr}_kaldi/cleantr.scp ark,t:${path}/${ftr}_kaldi/cleantr_cmvn.utt.ark
	done
	
done


