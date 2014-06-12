path=/home/li-bo/aurora2/data/noisemask/ftrs
curpath=`pwd`

for ftr in IBMFBankNoise TBMFBankNoise
do
	# training
	for data in cleantr multitr
	do
		cat ../fbank_e_d_a_kaldi/${data}.scp | awk '{print $1}' > tmp.scp
		copy-feats-from-text --dim=72 --data-suffix=txt --data-directory=ftrs/${ftr}/train tmp.scp ark,scp:${path}/${ftr}_kaldi/${data}.ark,${path}/${ftr}_kaldi/${data}.scp
	done
	
	# testing
	for data in testa testb testc
	do
		cat ../fbank_e_d_a_kaldi/${data}.scp | awk '{print $1}' > tmp.scp
		copy-feats-from-text --dim=72 --data-suffix=txt --data-directory=ftrs/${ftr}/${data} tmp.scp ark,scp:${path}/${ftr}_kaldi/${data}.ark,${path}/${ftr}_kaldi/${data}.scp
	done

	cd ${path}/${ftr}_kaldi
	compute-cmvn-stats scp:cleantr.scp ark,t:cleantr_cmvn.utt.ark
	compute-cmvn-stats scp:multitr.scp ark,t:multitr_cmvn.utt.ark
	compute-cmvn-stats scp:testa.scp ark,t:testa_cmvn.utt.ark
	compute-cmvn-stats scp:testb.scp ark,t:testb_cmvn.utt.ark
	compute-cmvn-stats scp:testc.scp ark,t:testc_cmvn.utt.ark
	cd ${curpath}

done


