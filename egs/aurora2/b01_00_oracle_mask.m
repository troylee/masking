clear all;
close all;

addpath(genpath('masking'));

resdir_ibm='ftrs/IBMOracle/';
resdir_tbm='ftrs/TBMOracle/';

data={'testa', 'testb', 'testc'}; % 

LC=-15;

%% Extract features
%% log Mel FBank domain masking
for i=1:length(data),
    
	disp(data(i));

    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
    for j=1:length(noisywav),
        
        %% IBMFBankDirect
        [noisyFBank, oriFBank, mask]=IBMFBankDirect_htk(char(noisywav(j)), char(cleanwav(j)), LC, 0);
        save([resdir_ibm char(tgt(j))], 'noisyFBank', '-ascii');
        
        %% TBMFBankDirect
        [noisyFBank, oriFBank, mask]=TBMFBankDirect_htk(char(noisywav(j)), char(cleanwav(j)), LC, 0);
        
        save([resdir_tbm char(tgt(j))], 'noisyFBank', '-ascii');
    end
end

    

