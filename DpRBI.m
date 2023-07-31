%% Description
% This MATLAB script computes the Dual-polarimteric Radar Built-up Index.
% It can be used for both VV-VH (Sentinel-1) and HH-HV (NISAR) pol modes.

%% Author Details
% Abhinav Verma, MRSLab, IIT Bombay
% Avik Bhattacharya, MRSLab, IIT Bombay

%% Related publication
% Abhinav Verma, Avik Bhattacharya, Subhadip Dey, Carlos López-Martínez, and Paolo Gamba 2023
% "Built-up Area Mapping Using Sentinel-1 SAR Data". ISPRS Journal of Photogrammetry and Remote Sensing.

%% List of inputs
% 1. Processed C2 (covraince matrix)
% 2. Slope (optional) -> This is to mask hilly and mountainous regions with high
% backscatter. It can be simply calculated in SNAP. Please note that if you
% do not have the slope image/data, turn off line no. 53 and 111 to 135
% However, we highly recommend to use slopes when applying DpRBI over
% hilly/mountainous regions.

%% List of outputs
% 1. Normalized decriptor g1n
% 2. Normalized decriptor g2n
% 3. Normalized decriptor g3n
% 4. DpRBI
% 5. Built-up Map

%% Reading the data

[filename, path] = uigetfile('*.*', 'Path selection');
path;
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp1 = fgets(f0);
tmp2 = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

%% C2 matrix

f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C22.bin'],'rb');
f5 = fopen([path 'slope.bin']); % Required for mountain mask

c11 = fread(f1,[ncols nrows],'float32') + ep;
c12 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
c21 = conj(c12);
c22 = fread(f4,[ncols nrows],'float32') + ep;
slope = fread(f5,[ncols nrows],'float32') + ep; %% if you turn off this, turn off line 111 to 135 (Mountain mask block)

fclose('all');

%% Intitialization
s0_convp = zeros(ncols,nrows);
s1_convp = zeros(ncols,nrows);
s2_convp = zeros(ncols,nrows);
s3_convp = zeros(ncols,nrows);
slope2 = zeros(ncols,nrows);
BA_flt = zeros(ncols,nrows);

%% for window processing

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing

for ii=startj:stopj
    for jj=starti:stopi

        %% C2 matrix

        c11c = mean2(c11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12c = mean2(c12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c21c = mean2(c21(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22c = mean2(c22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample

        C_2 = [c11c c12c; c21c c22c];

        % Stokes vector elements

        s0_convp(ii,jj) = (C_2(1,1) + C_2(2,2));
        s1_convp(ii,jj) = (C_2(1,1) - C_2(2,2));
        s2_convp(ii,jj) = 2*real(C_2(1,2));
        s3_convp(ii,jj) = 2*imag(C_2(1,2));

    end
    fprintf('Column: %d \n',ii);
end

%% Absolute Stokes values (only polarized componenets)

s1 = abs(s1_convp);
s2 = abs(s2_convp);
s3 = abs(s3_convp);

%% Mountain mask (comment off if mask is not required)

wsi_s = 9;
wsj_s = wsi_s; % Number of columns in the window

inci=fix(wsi_s/2); % Up & down movement margin from the central row
incj=fix(wsj_s/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti_s=fix(wsi_s/2)+1; % Starting row for window processing
startj_s=fix(wsj_s/2)+1; % Starting column for window processing

stopi_s= nrows-inci; % Stop row for window processing
stopj_s= ncols-incj; % Stop column for window processing

for ii=startj_s:stopj_s
    for jj=starti_s:stopi_s
        slope2(ii,jj) = mean2(slope(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
    end
end

slope_mask = slope2<=12;
slope_mask = double(slope_mask);

s1 = s1.*slope_mask;
s2 = s2.*slope_mask;
s3 = s3.*slope_mask;

%% Outlier Handling

% First polarized element
s1_vec = s1(:);
s1_5 = prctile(s1_vec,5);
s1_95 = prctile(s1_vec,95);
s1_cln = s1_vec(s1_vec > s1_5 & s1_vec < s1_95);
s1_max = max(s1_cln);
s1_min = min(s1_cln);
s1(s1 >=  s1_max) = s1_max;
s1(s1 <=  s1_min) = s1_min;

% Seceond polarized element
s2_vec = s2(:);
s2_5 = prctile(s2_vec,5);
s2_95 = prctile(s2_vec,95);
s2_cln = s2_vec(s2_vec > s2_5 & s2_vec < s2_95);
s2_max = max(s2_cln);
s2_min = min(s2_cln);
s2(s2 >=  s2_max) = s2_max;
s2(s2 <=  s2_min) = s2_min;

% Third polarized element
s3_vec = s3(:);
s3_5 = prctile(s3_vec,5);
s3_95 = prctile(s3_vec,95);
s3_cln = s3_vec(s3_vec > s3_5 & s3_vec < s3_95);
s3_max = max(s3_cln);
s3_min = min(s3_cln);
s3(s3 >=  s3_max) = s3_max;
s3(s3 <=  s3_min) = s3_min;

%% Normalized (Stokes) decriptors

s1_max2 = max(s1(:));
s1_norm = s1./s1_max2; %g1n

s2_max2 = max(s2(:));
s2_norm = s2./s2_max2; %g2n

s3_max2 = max(s3(:));
s3_norm = s3./(s3_max2); %g3n

%% Calculating DpRBI

DpRBI_11 = (sqrt(s1_norm.^2 + s2_norm.^2 + s3_norm.^2))./sqrt(3); %DpRBI

%% Built-up Area Map

dprbi_ba = DpRBI_11 >=0.58;

%% Post processing (BA Map Filtering using majority filtering)

wsi = 21;
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing

val2 = wsi*wsj;
for ii=startj:stopj
    for jj=starti:stopi

        ur = dprbi_ba(ii-inci:ii+inci,jj-incj:jj+incj);
        ur_1 = ur(:);
        val = sum(ur_1 == 1);

        vote = (val/val2)*100;
        if (vote >= 50)
            BA_flt(ii,jj) = 1;
        else
            BA_flt(ii,jj) = 0;
        end

    end
end

%% File Saving

f_name_100 = strcat(['g1n','.bin']);
fileandpath_100=strcat([path, f_name_100]);
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,s1_norm, 'float32');

f_name_200 = strcat(['g2n','.bin']);
fileandpath_200=strcat([path, f_name_200]);
fid_200 = fopen(fileandpath_200,'wb');
fwrite(fid_200,s2_norm, 'float32');

f_name_300 = strcat(['g3n','.bin']);
fileandpath_300=strcat([path, f_name_300]);
fid_300 = fopen(fileandpath_300,'wb');
fwrite(fid_300,s3_norm, 'float32');

f_name_400 = strcat(['DpRBI_index','.bin']);
fileandpath_400=strcat([path, f_name_400]);
fid_400 = fopen(fileandpath_400,'wb');
fwrite(fid_400,DpRBI_11, 'float32');

f_name_500 = strcat(['Built-up_Map','.bin']);
fileandpath_500=strcat([path, f_name_500]);
fid_500 = fopen(fileandpath_500,'wb');
fwrite(fid_500,BA_flt, 'float32');

fclose('all');



