clear all;
close all;

addpath(genpath('masking'));

respath='/home/li-bo/aurora2/exps/cleantr/x00_dbn/nnet/hid02d_dbn/decode';

resdir_ibm='ftrs/IBMFBankSpeech/';
resdir_tbm='ftrs/TBMFBankSpeech/';

data={'testa', 'testb', 'testc'};

IBM_LC=-6;
TBM_LC=0;

% load canonical clean speech spectral vectors
load('ftrs/cleanCanonicalVec.mat');

%% Extract features
for i=1:length(data),
    
	disp(data(i));

    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
	% load recognition results
	load([respath '/' char(data(i)) '.mat']);

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
end

    

