clear all;
close all;

addpath(genpath('masking'));

resdir_ibm='ftrs/IBMFBankSpeech/';
resdir_tbm='ftrs/TBMFBankSpeech/';

IBM_LC=-6;
TBM_LC=0;

% load canonical clean speech spectral vectors
load('ftrs/cleanCanonicalVec.mat');

%% Extract features
    
% load filenames
[noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_cleantr.scp'], '%s %s %s');
    
% load recognition results
load(['/home/li-bo/aurora2/exps/cleantr/a00_hmm/mono1b/pdf_align/train_pdf.mat']);

for j=1:length(noisywav),
	%% construct clean canonical speech based on recognition results
	cleanCanonicalSpec = ConstructSpeech(int32(pdfaligns{j}.states)+1, cleanCanonicalVec);
        
    %% IBMSpeechMask
    [maskedFBank, noisyPowerSpec, cleanPowerSpec, mask]=IBMSpeechMask(char(noisywav(j)), cleanCanonicalSpec, IBM_LC, 1);
        
    save([resdir_ibm char(tgt(j))], 'maskedFBank', '-ascii');
        
    %% TBMNoiseMask
    [maskedFBank, noisyPowerSpec, cleanPowerSpec, mask]=TBMSpeechMask(char(noisywav(j)), cleanCanonicalSpec, TBM_LC, 1);
        
    save([resdir_tbm char(tgt(j))], 'maskedFBank', '-ascii');
end

    

