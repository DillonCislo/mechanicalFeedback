mainDataDir = '/Users/idse/Dropbox/Ken/mechanical_feedback/';

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

i = i+1;
dataDir{i} = fullfile(mainDataDir, 'Ban_younger', '1d Ban clones');
bareFnames{i} = {   '20150221_Ban_JubGFP_1d_No.4' };
description{i} = 'Bantam 1.5d';
resolutions{i} = {0.110};
detectionChannel{i} = [1];
labels{i} = {'Jub', 'clone', 'arm'};

i = i+1;
dataDir{i} = fullfile(mainDataDir, 'Ban_younger', '1.5d Ban clones');
bareFnames{i} = {   '20150220_Ban_JubGFP_1.5d_No.2',...
                    '20150220_Ban_JubGFP_1.5d_No.3' };
description{i} = 'Bantam 1.5d';
resolutions{i} = {0.110, 0.105};
detectionChannel{i} = [1];
labels{i} = {'Jub', 'clone', 'arm'};

i = i+1;
dataDir{i} = fullfile(mainDataDir, 'Ban_younger','2d Ban clones');
bareFnames{i} = {   '20150221_Ban_JubGFP_2d_No.3',...
                    '20150226_Ban_JubGFP_2d_No.1'};
description{i} = 'Bantam 2d';
detectionChannel{i} = [1];
resolutions{i} = {0.110, 0.110};
labels{i} = {'Jub', 'clone', 'arm'};

%-------- bantam actin fullstack ---------------- 5

i = i+1;
dataDir{i} = fullfile(mainDataDir, 'Ban_actin_2.5d');
bareFnames{i} = {   '20150211; Ban_clones_Factin_DNA_No.1',...
                    '20150211; Ban_clones_Factin_DNA_No.4',...
                    '20150211; Ban_clones_Factin_DNA_No.5'};
description{i} = 'bantam actin';          
resolutions{i} = {0.134, 0.128, 0.121};
detectionChannel{i} = [2];
labels{i} = {'actin', 'arm', 'clone', 'DNA'};

%-------- Myc ---------------- 6

i = i+1;
dataDir{i} = fullfile(mainDataDir, 'Myc clones');
bareFnames{i} = {   '20150211_Myc_SqhGFP_No.5',...
                    %'20150211_Myc_SqhGFP_No.7',...,0.110
                    '20150211_Myc_SqhGFP_No.9'};
description{i} = 'Myc';   
resolutions{i} = { 0.110,0.100  };
detectionChannel{i} = [1];
labels{i} = {'Sqh', 'clone', 'arm'};

%-------- WT clones full stack ---------------- 7

i = i+1;
dataDir{i} = fullfile(mainDataDir, 'WTclones_SqhGFP');
bareFnames{i} = {   '20150720_WT_clones_SqhGFP_No.1' };
description{i} = 'WT clones';          
resolutions{i} = { 0.120 };
detectionChannel{i} = [1];
labels{i} = {'Sqh', 'DNA', 'clone', 'arm'};

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
resolutions{i} = { 0.120, 0.120, 0.120 };
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
bareFnames{i} = {   '20150603_Ban_E2F_RNAi_JubGFP_No.2',...
                    '20150603_Ban_E2F_RNAi_JubGFP_No.5'};
description{i} = 'Ban_E2F_RNAi';          
resolutions{i} = { 0.109, 0.109 };
detectionChannel{i} = [1];
labels{i} = {'Jub','clone', 'arm'};

