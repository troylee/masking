clear all;
close all;

addpath(genpath('masking'));

resdir_ibm='ftrs/IBMFBankNoise/';
resdir_tbm='ftrs/TBMFBankNoise/';

data={'train', 'testa', 'testb', 'testc'};

IBM_LC=0;
TBM_LC=0;

method='mcra';

%% Extract features
for i=1:length(data),
    
	disp(data(i));

    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
    for j=1:length(noisywav),
        
        %% IBMNoiseMask
        [maskedFBank, oriSpec, mask]=IBMNoiseMask(char(noisywav(j)), IBM_LC, 1, method);
        
        save([resdir_ibm char(tgt(j))], 'maskedFBank', '-ascii');
        
        %% TBMNoiseMask
        [maskedFBank, oriSpec, mask]=TBMNoiseMask(char(noisywav(j)), TBM_LC, 1, method);
        
        save([resdir_tbm char(tgt(j))], 'maskedFBank', '-ascii');
    end
end

    

