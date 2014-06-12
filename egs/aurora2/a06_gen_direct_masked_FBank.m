clear all;
close all;

addpath(genpath('masking'));

resdir_ibm='ftrs/IBMFBankDirect/';
resdir_tbm='ftrs/TBMFBankDirect/';

data={'train', 'testa', 'testb', 'testc'};

LC=0;

%% Extract features
for i=1:length(data),
    
	disp(data(i));

    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
    for j=1:length(noisywav),
        
        %% IBMFBankDirect
        [noisyFBank, oriFBank, mask]=IBMFBankDirect_htk(char(noisywav(j)), char(cleanwav(j)), LC, 1);
        
        save([resdir_ibm char(tgt(j))], 'noisyFBank', '-ascii');
        
        %% TBMFBankDirect
        [noisyFBank, oriFBank, mask]=TBMFBankDirect_htk(char(noisywav(j)), char(cleanwav(j)), LC, 1);
        
        save([resdir_tbm char(tgt(j))], 'noisyFBank', '-ascii');
    end
end

    

