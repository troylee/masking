respath=/home/li-bo/aurora2/exps/cleantr/a00_hmm/mono1b/eval_vts_model_n1e4

for data in testa testb testc
do
	cat ${respath}/${data}*.ali | ali-to-pdf ${respath}/../pdf_align/final.mdl ark:- ark,t:${respath}/${data}.pdf

	python pdf2mat.py ${data} ${respath}/${data}.pdf ${respath}/${data}.mat
done

