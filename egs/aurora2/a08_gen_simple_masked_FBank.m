clear all;
close all;

addpath(genpath('masking'));

resdir_ibm='ftrs/IBMFBankSimple/';
resdir_tbm='ftrs/TBMFBankSimple/';

data={'train', 'testa', 'testb', 'testc'};

IBM_LC=20;
TBM_LC=0;

%% Extract features
for i=1:length(data),
    
	disp(data(i));

    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
    for j=1:length(noisywav),
        
        %% IBMSimpleMask
        [maskedFBank, oriSpec, mask]=IBMSimpleMask(char(noisywav(j)), IBM_LC, 1);
        
        save([resdir_ibm char(tgt(j))], 'maskedFBank', '-ascii');
        
        %% TBMSimpleMask
        [maskedFBank, oriSpec, mask]=TBMSimpleMask(char(noisywav(j)), TBM_LC, 1);
        
        save([resdir_tbm char(tgt(j))], 'maskedFBank', '-ascii');
    end
end

    

