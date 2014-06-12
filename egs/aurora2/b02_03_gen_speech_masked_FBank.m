clear all;
close all;

addpath(genpath('masking'));

respath='/home/li-bo/aurora2/exps/cleantr/a00_hmm/mono1b/eval_vts_model_n1e4';

resdir_ibm='ftrs/IBMSpeech/';
resdir_tbm='ftrs/TBMSpeech/';

data={'testa', 'testb', 'testc'};

IBM_LC=-6;
TBM_LC=0;

% load canonical clean speech spectral vectors
load('ftrs/cleanPDFFBankVec.mat');

%% Extract features
for i=1:length(data),
    
	disp(data(i));

    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
	% load recognition results
	load([respath '/' char(data(i)) '.mat']);

    for j=1:length(noisywav),
		%% construct clean canonical speech based on recognition results
		cleanCanonicalSpec = ConstructSpeech(int32(pdfaligns{j}.states)+1, cleanPDFFBankVec);
        
        %% IBMSpeechMask
        %[maskedFBank, noisyPowerSpec, cleanPowerSpec, mask]=IBMSpeechMask(char(noisywav(j)), cleanCanonicalSpec, IBM_LC, 1);
        
        %save([resdir_ibm char(tgt(j))], 'maskedFBank', '-ascii');
        
        %% TBMNoiseMask
        [maskedFBank, oriFBank, mask]=TBMSpeechMask_FBank(char(noisywav(j)), cleanCanonicalSpec, TBM_LC, 1);
        
        save([resdir_tbm char(tgt(j))], 'maskedFBank', '-ascii');
    end
end

    

