mainDataDir = '/Users/idse/Dropbox/Ken/mechanical_feedback';


%-------- Bantam Sqh full stack ----------------

i = 1;
dataDir{i} = fullfile(mainDataDir, 'Ban_SqhGFP_2.5d');
bareFnames{i} = {  '20150610;Ban_clones_SqhGFP_No.3',...
                '20150610;Ban_clones_SqhGFP_No.4',...
                '20150610;Ban_clones_SqhGFP_No.5'};
resolutions{i} = { 0.110, 0.120, 0.120 };
description{i} = 'Bantam Squash full stack';
detectionChannel{i} = [3];
labels{i} = {'Sqh', 'DNA', 'clone', 'arm'};

%-------- Bantam Sqh younger half stack ----------------

dataDir{2} = fullfile(mainDataDir, 'Ban_younger', '1d Ban clones');
bareFnames{2} = {   '20150221_Ban_JubGFP_1d_No.4' };
description{2} = 'Bantam 1.5d';
resolutions{2} = {0.110};
detectionChannel{2} = [1];
labels{2} = {'Jub', 'clone', 'arm'};

dataDir{3} = fullfile(mainDataDir, 'Ban_younger', '1.5d Ban clones');
bareFnames{3} = {   '20150220_Ban_JubGFP_1.5d_No.2',...
                    '20150220_Ban_JubGFP_1.5d_No.3' };
description{3} = 'Bantam 1.5d';
resolutions{3} = {0.110, 0.105};
detectionChannel{3} = [1];
labels{3} = {'Jub', 'clone', 'arm'};

dataDir{4} = fullfile(mainDataDir, 'Ban_younger','2d Ban clones');
bareFnames{4} = {   '20150221_Ban_JubGFP_2d_No.3',...
                    '20150226_Ban_JubGFP_2d_No.1'};
description{4} = 'Bantam 2d';
detectionChannel{4} = [1];
resolutions{4} = {0.110, 0.110};
labels{4} = {'Jub', 'clone', 'arm'};

%-------- Myc ----------------

i = 5;
dataDir{i} = fullfile(mainDataDir, 'Myc clones');
bareFnames{i} = {   '20150211_Myc_SqhGFP_No.5',...
                    %'20150211_Myc_SqhGFP_No.7',...,0.110
                    '20150211_Myc_SqhGFP_No.9'};
description{i} = 'Myc';   
resolutions{i} = { 0.110,0.100  };
detectionChannel{i} = [1];
labels{i} = {'Sqh', 'clone', 'arm'};
            
%-------- WT clones full stack ----------------

dataDir{6} = fullfile(mainDataDir, 'WTclones_SqhGFP');
bareFnames{6} = {   '20150720_WT_clones_SqhGFP_No.1' };
description{6} = 'WT clones';          
resolutions{6} = { 0.120 };
detectionChannel{6} = [1];
labels{6} = {'Sqh', 'DNA', 'clone', 'arm'};

%-------- bantam actin fullstack ----------------

dataDir{7} = fullfile(mainDataDir, 'Ban_actin_2.5d');
bareFnames{7} = {   '20150211; Ban_clones_Factin_DNA_No.1',...
                    '20150211; Ban_clones_Factin_DNA_No.4',...
                    '20150211; Ban_clones_Factin_DNA_No.5'};
description{7} = 'bantam actin';          
resolutions{7} = {0.134, 0.128, 0.121};
detectionChannel{7} = [2];
labels{7} = {'actin', 'arm', 'clone', 'DNA'};

%-------- WT in minute full stack ----------------

i = 8;
dataDir{i} = fullfile(mainDataDir, 'WTinMinute_ZipGFP_2.5d');
bareFnames{i} = {   '20150618_WT_in_Min_ZipGFP_No.2',...
                    '20150618_WT_in_Min_ZipGFP_No.3'};
description{i} = 'WT in minute';          
resolutions{i} = { 0.120, 0.120 };
detectionChannel{i} = [1];
labels{i} = {'Sqh', 'clone', 'DNA', 'arm'};

%-------- ROK rnai ---------------- 9
i = i + 1;
dataDir{i} = fullfile(mainDataDir, 'RhoKRNAi_SqhGFP_2.5d');
bareFnames{i} = {   '20150720_RhoK_RNAi_SqhGFP_No.1',...
                    '20150720_RhoK_RNAi_SqhGFP_No.3',...
                    '20150720_RhoK_RNAi_SqhGFP_No.4'};
description{i} = 'Rok RNAi';          
resolutions{i} = { 0.120,0.120,0.120 };
detectionChannel{i} = [1];
labels{i} = {'Sqh', 'DNA', 'clone', 'arm'};

%-------- pure WT ---------------- 10

i = i + 1;
dataDir{i} = fullfile(mainDataDir, 'WT_SqhGFP');
bareFnames{i} = {   '20150717_WT_SqhGFP_No.1',...
                    '20150717_WT_SqhGFP_No.2'};
description{i} = 'WT';          
resolutions{i} = { 0.120,0.120 };
detectionChannel{i} = [1];
labels{i} = {'Sqh', 'DNA', 'arm'};

%-------- Ban E2F RNAi ---------------- 11

i = i + 1;
dataDir{i} = fullfile(mainDataDir, 'Ban_E2F_RNAi_Jub_2.5d');
bareFnames{i} = {   '20150603_Ban_E2F_RNAi_JubGFP_No.2'};
description{i} = 'Ban_E2F_RNAi';          
resolutions{i} = { 0.109 };
detectionChannel{i} = [1];
labels{i} = {'Jub','clone', 'arm'};


%--------- Ban JubGFP 2.5d -----------12
i = i + 1;
dataDir{i} = fullfile(mainDataDir, '20150603_Ban_JubGFP_2.5d');
bareFnames{i} = {   '20150603_Ban_JubGFP_2.5d_No.2'};
description{i} = 'Ban_JubGFP';          
resolutions{i} = { 0.109 };
detectionChannel{i} = [1];
labels{i} = {'Jub','clone', 'arm'};
