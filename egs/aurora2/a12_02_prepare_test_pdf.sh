respath=/home/li-bo/aurora2/exps/cleantr/x00_dbn/nnet/hid02d_dbn/decode

for data in testa testb testc
do
	cat ${respath}/${data}*.ali | ali-to-pdf ${respath}/../../../lib/final.mdl ark:- ark,t:${respath}/${data}.pdf

	python pdf2mat.py ${data} ${respath}/${data}.pdf ${respath}/${data}.mat
done

