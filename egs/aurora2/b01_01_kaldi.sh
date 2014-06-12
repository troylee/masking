path=/home/li-bo/aurora2/data/noisemask/ftrs
curpath=`pwd`

train(){
	# training
    for data in cleantr multitr
    do
        cat ../fbank_e_d_a_kaldi/${data}.scp | awk '{print $1}' > tmp.scp
        copy-feats-from-text --dim=24 --data-suffix=txt --data-directory=ftrs/${ftr}/train tmp.scp ark,scp:${path}/${ftr}_kaldi/${data}.ark,${path}/${ftr}_kaldi/${data}.scp

		compute-cmvn-stats scp:${path}/${ftr}_kaldi/${data}.scp ark,t:${path}/${ftr}_kaldi/${data}_cmvn.utt.ark
    done
}

test(){
	# testing
    for data in testa testb testc
    do
        cat ../fbank_e_d_a_kaldi/${data}.scp | awk '{print $1}' > tmp.scp
        copy-feats-from-text --dim=24 --data-suffix=txt --data-directory=ftrs/${ftr}/${data} tmp.scp ark,scp:${path}/${ftr}_kaldi/${data}.ark,${path}/${ftr}_kaldi/${data}.scp

		compute-cmvn-stats scp:${path}/${ftr}_kaldi/${data}.scp ark,t:${path}/${ftr}_kaldi/${data}_cmvn.utt.ark
    done
}

for ftr in IBMOracle TBMOracle
do

	#train
	test
done


