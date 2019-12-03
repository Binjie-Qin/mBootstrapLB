%%  Generate a random sample index sequence with bootstrap resample method

[~, bootsamp] = bootstrp(nboot, bootfun, dl);

%   INPUT ARGUMENTS:
%       1) nboot: Generate nboot bootstrap data samples. 
%                  (In this project, nboot = 300);
%       2) bootfun: Compute statistics on each sample using bootfun.
%                  (In this project, bootfun = @mean);
%       3) dl: The data (scalars, column vectors, or matrices) used to
%               create inputs to bootfun.
%                   (In this project, dl = 0:8)
%
%
%   SO, in this project, the bootstrap resample operator is:
%
%       [~, bootsamp] = bootstrp(300, @mean, 0:8);
%
%
%%  Compute Poisson regression using bootstrap method

A = mBootstrapLB(bootsamp, pp_data, p_data);

%   INPUT ARGUMENTS:
%       1) bootsamp: The random sample index sequence using bootstrap
%                       resample method (The result of the above equation);
%       2) pp_data: The origin image data to be processed;
%       3) p_data: The logarithmic form of above the image data.
%
%   OUTPUT ARGUMENTS:
%       A: The logarithmic form of denoising image data
%
%
%%  BPR step:

p_data = data;
pp_data = log(p_data + 0.5);
[~, bootsamp] = bootstrp(300, @mean, 0:8);
B = mBootstrapLB(bootsamp, pp_data, p_data);
A = exp(B) - 0.5;

%   PARAMETER DECLARATION:
%       1) data: The origin image data;
%       2) p_data: The origin image data;
%       3) pp_data: The logarithmic form of the origin image data;
%       4) bootsamp: The random sample index sequence;
%       5) B: The logarithmic form of the denoising image data;
%       6) A: The denoising iamge data after processing.
%
%
%%  NPR algorithm:

B = NPR(X,Y,h);
%   INPUT ARGUMENTS:
%       1) X: The X coordinate sequence corresponding to the spectrum
%               domain (i.e., if there are 1024 channels in the detector,
%               X must be the incresing sequence such as 1:1024);
%               (In this project, X = 1:1024 for detecting the metal data);
%       2) Y: The spectrum data which need to calculate the baseline;
%       3) h: The range of the local regression;
%               (In this project, h = 250).
%
%   OUTPUT ARGUMENTS:
%       B: the baseline of the spectrum data.
%
%
%%  Remove the baseline algorithm (RR step):

% Example:

X = 1:1024;
X = X';
D = data(:,x);
baseline = NPR(X, D, 250);
result = D - baseline;

%   PARAMETER DECLARATION:
%       1) X: The X coordinate sequence corresponding to the spectrum
%               domain (i.e., if there are 1024 channels in the detector,
%               X must be the incresing sequence such as 1:1024);
%       2) data: The origin image data;
%       3) D: The spectrum data which need to calculate the baseline;
%       4) x: x present a certain point which has high intensity;
%       5) baseline: the baseline of the spectrum data.
%       6) result: The actual spectrum after removing the baseline.
%
%