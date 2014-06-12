%
% accumulate clean spectral features to estimate the clean canonical
% vectors.
%

clear all;
close all;

addpath(genpath('masking'));

numpdfs=181;

specdim=129;

numfiles=8440;

% load the pdf alignments for all the clean train data, pdfaligns
load /home/li-bo/aurora2/exps/cleantr/a00_hmm/mono1b/pdf_align/train_pdf.mat;

% load the scp file list
[wavfiles, cleanwavs, tgts]=textread('/home/li-bo/aurora2/data/lists/wav2mat_train_all.scp', '%s %s %s');

cleanCanonicalVec=zeros(numpdfs, specdim);
vecCount=zeros(numpdfs,1);

fidx=1;
for i=1:length(wavfiles),
    if fidx<=numfiles && strcmp(char(wavfiles(i)), ['train/' pdfaligns{fidx}.name '.txt'])==0,
        
        disp([num2str(fidx) ' --- ' char(wavfiles(i))]);
        
        specftr=Spectrum_htk(char(wavfiles(i)), 0);

		% in Kaldi the alignment pdf id starts from 0
		states=int32(pdfaligns{fidx}.states) + 1;
        
        if size(specftr,1) ~= length(states),
            disp('Unmatched number of frames!');
            return;
        end
        if size(specftr,2) ~= specdim,
            disp('Spectral feature dimension error!');
            return;
        end
        
        for k=1:length(states),
            pdfid=states(k);
            
            vecCount(pdfid)=vecCount(pdfid)+1;
            cleanCanonicalVec(pdfid,:)=cleanCanonicalVec(pdfid,:)+specftr(k,:);
        end
        fidx=fidx+1;
    end
end

% compute the mean
for i=1:numpdfs,
    cleanCanonicalVec(i,:)=cleanCanonicalVec(i,:)/vecCount(i);
end

% save the results
save('ftrs/cleanCanonicalVec.mat', 'cleanCanonicalVec');


