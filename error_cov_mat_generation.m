%%% generate errors covariance matrix for weight determination of
%%% observations in spowtd PEST calibration for CO PEATLCSM parameters

% read in the file with the observations for Ekolongouma and Itanga
delimiterIn = ' '; %define seperator
headerlinesIn = 1; %define header
filename_itanga = '/data/leuven/324/vsc32460/AC/spowtd/curves_calibration_poros_shift_itanga_cov.res';
filename_ekolongouma = '/data/leuven/324/vsc32460/AC/spowtd/curves_calibration_poros_shift_cov.res';
Itanga = importdata(filename_itanga,delimiterIn,headerlinesIn); %import the data
Ekolongouma = importdata(filename_ekolongouma,delimiterIn,headerlinesIn); %import the data

Itanga_obs_stor = Itanga.data(1:841,3);
Itanga_obs_time = Itanga.data(842:end,3);

Ekolongouma_obs_stor = Ekolongouma.data(1:689,3);
Ekolongouma_obs_time = Ekolongouma.data(690:end,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    ADJUST THE FOLLOWING TWO VARIABLES                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% either 'time' or 'storage'
observation = "time";
% either 'ekolongouma' or 'itanga'
site = "ekolongouma";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(site,'ekolongouma')
    N_obs_stor = length(Ekolongouma_obs_stor);
    N_obs_time = length(Ekolongouma_obs_time);
    pd = fitdist(Ekolongouma_obs_time,'Normal');
    sigma_obs_time = pd.sigma;
    pd = fitdist(Ekolongouma_obs_stor,'Normal');
    sigma_obs_stor =pd.sigma;
elseif strcmp(site,'itanga')
    N_obs_stor = length(Itanga_obs_stor);
    N_obs_time = length(Itanga_obs_time);
    pd = fitdist(Itanga_obs_time,'Normal');
    sigma_obs_time = pd.sigma;
    pd = fitdist(Itanga_obs_stor,'Normal');
    sigma_obs_stor =pd.sigma;
end
    
% in mm storage obs and in days timeobs
sigma_obs_stor = 8.0;    %7.2     %choose sd of observation error for which the cov matrix is calculated (storageobs/timeobs!)
sigma_obs_time = 2.0 ;       %2.5

%error corr length (number of samples/observations that are autocorrelated)
rho_stor = 50;
rho_time = 50;

if strcmp(observation,'time')
N_obs=N_obs_time;
sigma_obs = sigma_obs_time;
rho = rho_time;
elseif strcmp(observation,'storage')
N_obs=N_obs_stor;
sigma_obs = sigma_obs_stor;
rho = rho_stor;
end

M=nan+zeros(N_obs,N_obs);

for i = 1:N_obs
    for j = i:N_obs
        x=j-i;
        r = exp(-x/rho);
        M(i,j)=r*sigma_obs.^2;
        M(j,i)=M(i,j);
    end
end

%M
outpath = '/data/leuven/324/vsc32460/AC/spowtd/';
fullfile_name = fullfile(outpath, string(site) +'_' +string(observation) +'obs_errcovD.dat');
save(fullfile_name, 'M','-ascii','-double')
%saveas(M,fullfile_name,'Delimiter','space');
exit

