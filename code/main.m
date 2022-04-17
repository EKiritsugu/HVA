%% load data

misc.dataFolder = '../data'; % This place is the requirement of Code Ocean.

misc.dataName1 = 'dev1_female4_liverec_130ms_5cm_sim_1.wav';
misc.dataName2 = 'dev1_female4_liverec_130ms_5cm_sim_4.wav';

if ~exist([misc.dataFolder '/' misc.dataName1],'file')
    if strcmp('OK!',questdlg(['Can I download 200+ wav files to "' misc.dataFolder '"?'],'Downloading dataset','OK!','NO!','NO!'))
        disp 'Downloading data...'
        unzip('http://www.irisa.fr/metiss/SiSEC10/underdetermined/dev1.zip',misc.dataFolder)
        disp 'Audio files were downloaded.'
    else
        error(['Specified data files do not exist in "' misc.dataFolder '".'])
    end
end

[signal1,~       ] = audioread([misc.dataFolder '/' misc.dataName1]);
[signal2,param.fs] = audioread([misc.dataFolder '/' misc.dataName2]);


%% set parameters

param.STFTwinSec = 0.128;
param.STFTwinShift = 2;

param.STFTwindowLength = round(param.STFTwinSec*param.fs);
param.STFTshiftSize = param.STFTwindowLength/param.STFTwinShift;

param.sourceNum = 2;
param.refMicIndex = 1;

param.mu1 = 1;
param.mu2 = 1/param.mu1;
param.alpha = 1;
param.iterNum = 100;

param.lambda = 0.08;
param.kappa = 3;


%% mix and separate

mixture = signal1 + signal2;

tic
separated = HVA(mixture,param);
toc


%% plot results

h = figure;

subplot(5,1,1)
plot(signal1(:,1))
grid on
axis tight
xticklabels([])
title 'Signal 1'

subplot(5,1,2)
plot(signal2(:,1))
grid on
axis tight
xticklabels([])
title 'Signal 2'

subplot(5,1,3)
plot(mixture)
grid on
axis tight
xticklabels([])
title 'Mixture'

subplot(5,1,4)
plot(separated(:,1))
grid on
axis tight
xticklabels([])
title 'Separated signal'

subplot(5,1,5)
plot(separated(:,2))
grid on
axis tight
xticklabels([])
title 'Separated signal'

saveas(h,'../results/result.png')