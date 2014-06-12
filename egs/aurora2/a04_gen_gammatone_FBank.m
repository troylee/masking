clear all;
close all;

addpath(genpath('masking'));

resdir='ftrs/gammafbk/';

data={'train', 'testa', 'testb', 'testc'}; % 'train'

useDynamic=1;

%% Extract features
for i=1:length(data),
	
	disp(data(i));
    
    % load filenames
    [noisywav, cleanwav, tgt]=textread(['../lists/wav2mat_' char(data(i)) '_all.scp'], '%s %s %s');
    
    for j=1:length(noisywav),
        
        fbks = GammatoneFBank(char(noisywav(j)), useDynamic);
		
        save([resdir char(tgt(j))], 'fbks', '-ascii');

    end
end

    

