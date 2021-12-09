% This is a Matlab script to simulate sampling and quantization errors. To represent 
% the analog signal, a large sampling rate (4000 Hz) is used to represent continuous time. 
% Also, the numerical error is too tiny based on which we can
% assume amplitude is continuous too.

% Note 1: A uniform spectral density is assumed for the original analog signal
% Note 2: Amplitude is less than the clipping level (no satuartion error)
% Note 3: "Sampling" error is the error between S_d and S_c
% Note 4: "Sampling+Quantization" error is the error between S_D and S_c

% convention:
%       _c: analog continuous-time signal
%       _d: Discrete-time signal (continuous-amplitude)
%       _D: Digital signal (discrete-time and discrete-amplitude)    

% Written by S. Farid Ghahari (ghahari@gmail.com), 3/7/2021
% Revised by Wenjie Liao (liaowj17@mails.tsinghua.edu.cn), 5/4/2021

clear
clc
close all
rng(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%% read and resave acceleration data with noise  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read GM names
GM_root_Path = '.\0_PEER_Database\';  % raw data path
NGA_Flatfile_dir = [GM_root_Path,'NGA_Flatfile.xls']; % the Excel file recorded ground motion information
[rawDataNum, rawDataStr] = xlsread(NGA_Flatfile_dir); % read all info in NGA_Flatfile.xls
[num_GMs,c] = size(rawDataNum); % num_GMs is the number of GMs
sensor = 'Episensor'; % name of sensor [Phone, Phidgets, Episensor]

%% create root dir
PEER_RE_root = ['.\',sensor,'_PEER_Re\']; % ADC resample data
if ~exist(PEER_RE_root,'dir')
        mkdir(PEER_RE_root);
end
PEER_noise_root = ['.\',sensor,'_PEER_Noise\']; % ADC plus noise data
if ~exist(PEER_noise_root,'dir')
        mkdir(PEER_noise_root);
end
Analysis_root = ['.\',sensor,'_PEER_Analysis\']; % ADC plus noise data
if ~exist(Analysis_root,'dir')
        mkdir(Analysis_root);
end

%% set number groups
SNRs = zeros(num_GMs,1);
RMSEs = zeros(num_GMs,1);

%% ADC for PEER ground motions
for num_GM = 1:4000
    NGA_GM_num = ['NGA',rawDataStr{(num_GM+1),2}]; % get GM number
    NGA_GM_Mg = rawDataNum((num_GM),10); % get GM magnitude
    str_NGA_GM_Mg = num2str(NGA_GM_Mg); %cell2mat(regexp(num2str(double(NGA_GM_Mg)),'\d+.\d{1}','match'));
    NGA_GM_EpiD = rawDataNum((num_GM),48); % get GM epicenter distance
    str_NGA_GM_EpiD = cell2mat(regexp(num2str(double(NGA_GM_EpiD)),'\d+.\d{1}','match'));
    NGA_GM_SaT1 = rawDataNum((num_GM),196); % get GM SaT1 data
    %str_NGA_GM_SaT1 = regexp(num2str(NGA_GM_SaT1),'\d+.\d{3}','match');
    
	% enter file
    GM_current_Path = [GM_root_Path,'data\',NGA_GM_num,'\'];

    GM_RE_Path = [PEER_RE_root,'Mg',str_NGA_GM_Mg,'_EpiD',str_NGA_GM_EpiD,'_',NGA_GM_num,'\']; % create output path
    if ~exist(GM_RE_Path,'dir')
        mkdir(GM_RE_Path)
    end
    GM_noise_Path = [PEER_noise_root,'Mg',str_NGA_GM_Mg,'_EpiD',str_NGA_GM_EpiD,'_',NGA_GM_num,'\']; % create output path
    if ~exist(GM_noise_Path,'dir')
        mkdir(GM_noise_Path)
    end
    GM_Analysis_path = [Analysis_root,'Mg',str_NGA_GM_Mg,'_EpiD',str_NGA_GM_EpiD,'_',NGA_GM_num,'\']; % create output path
    if ~exist(GM_Analysis_path,'dir')
        mkdir(GM_Analysis_path)
    end
    
    AT_files = dir(fullfile(GM_current_Path,'*.AT2'));  % show files with .AT2
    AT_filenames = {AT_files.name}';            % get file names and reshape to 1 column 

    for i = 1:length(AT_filenames)
         % read data
        file_dir = strcat(GM_current_Path, AT_filenames(i));
        [NPTS,dt,~,~] = textread(file_dir{1,1},'%d %f %s %s',1,'headerlines', 3); % read basic info
        [col1,col2,col3,col4,col5] = textread(file_dir{1,1},'%f %f %f %f %f','headerlines',4);
        temp_Acc = zeros(length(NPTS),1);
        
        k = 1;
        for j = 1:length(col5)
            temp_Acc(k,1) = col1(j)*9.8; %acceleration x series,¡ÁPGA_scale, unit£ºm/s/s
            temp_Acc(k+1,1) = col2(j)*9.8; %acceleration x series,¡ÁPGA_scale, unit£ºm/s/s
            temp_Acc(k+2,1) = col3(j)*9.8; %acceleration x series,¡ÁPGA_scale, unit£ºm/s/s
            temp_Acc(k+3,1) = col4(j)*9.8; %acceleration x series,¡ÁPGA_scale, unit£ºm/s/s
            temp_Acc(k+4,1) = col5(j)*9.8; %acceleration x series,¡ÁPGA_scale, unit£ºm/s/s
            k = k+5;
        end
        
        NPTS = k-1; % update NPTS
        Acc = temp_Acc(1:NPTS,1);
        temp_time = (0:dt:NPTS*dt)';
        A_times = temp_time(1:NPTS,1);
        A_acc_xs = zeros(length(NPTS),1);
        A_acc_noise = zeros(length(NPTS),1);
        
        %% ADC process
        % for acceleration x
        fmax_c = 1/dt;
        T = NPTS*dt;
        [Time_d,A_acc_noise,A_acc,length_Acc,SNR] = ADC_with_noise_0519( A_times,Acc,fmax_c,T );
        A_acc_xs_noise_out = zeros(length_Acc,2);
        A_acc_xs_noise_out(:,1) = Time_d(1:length_Acc,1);
        A_acc_xs_noise_out(:,2) = A_acc_noise(1:length_Acc,1);

        A_acc_xs_re = zeros(length_Acc,2);
        A_acc_xs_re(:,1) = Time_d(1:length_Acc,1);
        A_acc_xs_re(:,2) = A_acc(1:length_Acc,1);

        SNRs(i,1) = SNR; % record SNR

        % resave signal data
        file_x_dir_re = strcat(GM_RE_Path, AT_filenames(i));
        file_x_dir_re = strrep(file_x_dir_re{1,1},'.AT2','.txt');
        % save as txt file
        fileID_re = fopen(file_x_dir_re,'w'); % file id
        num_data = length(A_acc_xs_re);
        fprintf (fileID_re,'%d\n',num_data);
        for k=1:num_data
            fprintf(fileID_re,'%.5f \t %.10f\n',A_acc_xs_re(k,:));
        end
        fclose(fileID_re);
        % save as CSV file
        csvfile_x_dir_re = strrep(file_x_dir_re,'.txt','.csv');
        csvwrite(csvfile_x_dir_re,A_acc_xs_re);% writes matrix M into FILENAME as comma-separated values.

        % save signal+noise data
        file_x_dir = strcat(GM_noise_Path, AT_filenames(i));
        file_x_dir = strrep(file_x_dir{1,1},'.AT2','.txt');
        % save as txt file
        fileID = fopen(file_x_dir,'w'); % file id
        num_data = length(A_acc_xs_noise_out);
        fprintf (fileID,'%d\n',num_data);
        for k=1:num_data
            fprintf(fileID,'%.5f \t %.10f\n',A_acc_xs_noise_out(k,:));
        end
        fclose(fileID);
        % save as CSV file
        csvfile_x_dir = strrep(file_x_dir,'.txt','.csv');
        csvwrite(csvfile_x_dir,A_acc_xs_noise_out);% writes matrix M into FILENAME as comma-separated values.

        %% spectrum analysis
        % no noise
        [Ts_x,Disp_x,Vel_x,AbsAcce_x] = Spectrum_0519(A_acc_xs_re);
        AbsAcce_x_out = zeros(length(AbsAcce_x),2);
        AbsAcce_x_out(:,1) = Ts_x';
        AbsAcce_x_out(:,2) = AbsAcce_x';
        % with noise 
        [Ts_x,Disp_x_noise,Vel_x_noise,AbsAcce_x_noise] = Spectrum_0519(A_acc_xs_noise_out);
        AbsAcce_x_noise_out = zeros(length(AbsAcce_x_noise),2);
        AbsAcce_x_noise_out(:,1) = Ts_x';
        AbsAcce_x_noise_out(:,2) = AbsAcce_x_noise';

        % save spectrum data
        % no noise
        if ~exist([GM_Analysis_path,'ADC_spect_',sensor,'\'],'dir')
            mkdir([GM_Analysis_path,'ADC_spect_',sensor,'\'])
        end
        file_spect_x_dir = strcat([GM_Analysis_path,'ADC_spect_',sensor,'\'], AT_filenames(i));
        save(file_spect_x_dir{1,1},'AbsAcce_x_out','-ascii');
        % with noise
        if ~exist([GM_Analysis_path,'ADC_spect_noise_',sensor,'\'],'dir')
            mkdir([GM_Analysis_path,'ADC_spect_noise_',sensor,'\'])
        end
        file_spect_x_dir = strcat([GM_Analysis_path,'ADC_spect_noise_',sensor,'\'], AT_filenames(i));
        save(file_spect_x_dir{1,1},'AbsAcce_x_noise_out','-ascii');

        %% spectrum RMSE calculate 
        error_AbsAcce_x = AbsAcce_x - AbsAcce_x_noise;
        RMSE_AbsAcce_x = sqrt(mean(sum(error_AbsAcce_x.^2)));
        RMSEs(i,1) = RMSE_AbsAcce_x;

    end
end

% save SNR and RMSE data
file_SNR_x_dir = [Analysis_root,'SNR_',sensor,'.txt'];
save(file_SNR_x_dir,'SNRs','-ascii');
file_RMSE_x_dir = [Analysis_root,'RMSE_',sensor,'.txt'];
save(file_RMSE_x_dir,'RMSEs','-ascii');
