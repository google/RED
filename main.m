% Copyright 2017 Google Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     https://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Demonstration of the image restoration experiments conducted in
% Y. Romano, M. Elad, and P. Milanfar, "The Little Engine that Could: 
% Regularization by Denoising (RED)", submitted to SIAM Journal on Imaging
% Sciences, 2016. https://arxiv.org/abs/1611.02862
%
% This example reads a ground-truth image, degrades the image by 
% first blurring or downscaling it, followed by an addition of random white
% Gaussian noise. Then it calls to RED in order to restore the image. 
% This example compares the input and output PSNR, shows and saves the 
% results. The suggested image-adaptive Laplacian-regularization functional 
% is minimized using the Fixed-Point, ADMM, and Steepest Descent methods. 
% Please refer to the paper for more details.
%
% The following are the degradation models that this example handles:
% 'UniformBlur'  - 9X9 uniform psf with noise-level equal to sqrt(2)
% 'GaussianBlur' - 25X25 Gaussian psf with std 1.6 and noise-level
%                  equal to sqrt(2)
% 'Downscale'    - 7X7 Gaussian psf with std 1.6 and noise-level equal to 5
%
% The denoising engine is TNRD: Yunjin Chen, and Thomas Pock, "Trainable 
% Nonlinear Reaction Diffusion: A Flexible Framework for Fast and Effective
% Image Restoration", IEEE TPAMI 2016. The code is available in
% http://www.icg.tugraz.at/Members/Chenyunjin/about-yunjin-chen
% Note: Enable parallel pool to reduce runtime.
%
% The degradation process is similar to the one suggested in NCSR paper:
% Weisheng Dong, Lei Zhang, Guangming Shi, and Xin Li "Nonlocally 
% Centralized Sparse Representation for Image Restoration", IEEE-TIP, 2013.
% The code is available in http://www4.comp.polyu.edu.hk/~cslzhang/NCSR.htm
%

clc;
clear;
close all;

% configure the path
% denoising functions
addpath(genpath('./tnrd_denoising/'));
% SD, FP, and ADMM methods
addpath(genpath('./minimizers/'));
% contains the default params
addpath(genpath('./parameters/'));
% contains basic functions
addpath(genpath('./helper_functions/'));
% test images for the debluring and super resolution problems, 
% taken from NCSR software package
addpath(genpath('./test_images/'));

% set light_mode = true to run the code in a sub optimal but faster mode
% set light_mode = false to obtain the results reported in the RED paper
light_mode = false;

if light_mode
    fprintf('Running in light mode. ');
    fprintf('Turn off to obatain the results reported in RED paper.\n');
else
    fprintf('Light mode option is off. ');
    fprintf('Reproducing the result in RED paper.\n');
end

%% read the original image

file_name = 'starfish.tif';

fprintf('Reading %s image...', file_name);
orig_im = imread(['./test_images/' file_name]);
orig_im = double(orig_im);

fprintf(' Done.\n');


%% define the degradation model

% choose the secenrio: 'UniformBlur', 'GaussianBlur', or 'Downscale'
degradation_model = 'UniformBlur';

fprintf('Test case: %s degradation model.\n', degradation_model);

switch degradation_model
    case 'UniformBlur'
        % noise level
        input_sigma = sqrt(2);
        % filter size
        psf_sz = 9;
        % create uniform filter
        psf = fspecial('average', psf_sz);
        % use fft to solve a system of linear equations in closed form
        use_fft = true;
        % create a function-handle to blur the image
        ForwardFunc = ...
            @(in_im) imfilter(in_im,psf,'conv','same','circular');
        % the psf is symmetric, i.e., the ForwardFunc and BackwardFunc
        % are the same
        BackwardFunc = ForwardFunc;
        % special initialization (e.g. the output of other method)
        % set to identity mapping
        InitEstFunc = @(in_im) in_im;
        
    case 'GaussianBlur'
        % noise level
        input_sigma = sqrt(2);
        % filter size
        psf_sz = 25;
        % std of the Gaussian filter
        gaussian_std = 1.6;
        % create gaussian filter
        psf = fspecial('gaussian', psf_sz, gaussian_std);
        % use fft to solve a system of linear equations in closed form
        use_fft = true;
        % create a function handle to blur the image
        ForwardFunc = ...
            @(in_im) imfilter(in_im,psf,'conv','same','circular');
        % the psf is symmetric, i.e., the ForwardFunc and BackwardFunc
        % are the same
        BackwardFunc = ForwardFunc;
        % special initialization (e.g. the output of other method)
        % set to identity mapping
        InitEstFunc = @(in_im) in_im;
        
    case 'Downscale'
        % noise level
        input_sigma = 5;
        % filter size
        psf_sz = 7;
        % std of the Gaussian filter
        gaussian_std = 1.6;
        % create gaussian filter
        psf = fspecial('gaussian', psf_sz, gaussian_std);
        % scaling factor
        scale = 3;
        
        % compute the size of the low-res image
        lr_im_sz = [ceil(size(orig_im,1)/scale),...
                    ceil(size(orig_im,2)/scale)];        
        % create the degradation operator
        H = CreateBlurAndDecimationOperator(scale,lr_im_sz,psf);
        % downscale
        ForwardFunc = @(in_im) reshape(H*in_im(:),lr_im_sz);        
        % upscale
        BackwardFunc = @(in_im) reshape(H'*in_im(:),scale*lr_im_sz);
        % special initialization (e.g. the output of other method)
        % use bicubic upscaler
        InitEstFunc = @(in_im) imresize(in_im,scale,'bicubic');
        
    otherwise
        error('Degradation model is not defined');
end


%% degrade the original image

switch degradation_model
    case {'UniformBlur', 'GaussianBlur'}
        fprintf('Blurring...');
        % blur each channel using the ForwardFunc
        input_im = zeros( size(orig_im) );
        for ch_id = 1:size(orig_im,3)
            input_im(:,:,ch_id) = ForwardFunc(orig_im(:,:,ch_id));
        end
        % use 'seed' = 0 to be consistent with the experiments in NCSR
        randn('seed', 0);

    case 'Downscale'
        fprintf('Downscaling...');
        % blur the image, similar to the degradation process of NCSR
        input_im = Blur(orig_im, psf);
        % decimate
        input_im = input_im(1:scale:end,1:scale:end,:);
        % use 'state' = 0 to be consistent with the experiments in NCSR
        randn('state', 0);

    otherwise
        error('Degradation model is not defined');
end

% add noise
fprintf(' Adding noise...');
input_im = input_im + input_sigma*randn(size(input_im));

% convert to YCbCr color space if needed
input_luma_im = PrepareImage(input_im);
orig_luma_im = PrepareImage(orig_im);

if strcmp(degradation_model,'Downscale')
    % upscale using bicubic
    input_im = imresize(input_im,scale,'bicubic');
    input_im = input_im(1:size(orig_im,1), 1:size(orig_im,2), :); 
end
fprintf(' Done.\n');
psnr_input = ComputePSNR(orig_im, input_im);


%% minimize the Laplacian regularization functional via Fixed Point

fprintf('Restoring using RED: Fixed-Point method\n');

switch degradation_model
    case 'UniformBlur'
        params_fp = GetUniformDeblurFPParams(light_mode, psf, use_fft);
    case 'GaussianBlur'
        params_fp = GetGaussianDeblurFPParams(light_mode, psf, use_fft);
    case 'Downscale'
        assert(exist('use_fft','var') == 0);
        params_fp = GetSuperResFPParams(light_mode);
    otherwise
        error('Degradation model is not defined');
end

[est_fp_im, psnr_fp] = RunFP(input_luma_im,...
                             ForwardFunc,...
                             BackwardFunc,...
                             InitEstFunc,...
                             input_sigma,...
                             params_fp,...
                             orig_luma_im);
out_fp_im = MergeChannels(input_im,est_fp_im);

fprintf('Done.\n');


%% minimize the Laplacian regularization functional via ADMM

fprintf('Restoring using RED: ADMM method\n');

switch degradation_model
    case 'UniformBlur'
        params_admm = GetUniformDeblurADMMParams(light_mode, psf, use_fft);
    case 'GaussianBlur'
        params_admm = GetGaussianDeblurADMMParams(light_mode, psf, use_fft);
    case 'Downscale'
        assert(exist('use_fft','var') == 0);
        params_admm = GetSuperResADMMParams(light_mode);        
    otherwise
        error('Degradation model is not defined');
end

[est_admm_im, psnr_admm] = RunADMM(input_luma_im,...
                                   ForwardFunc,...
                                   BackwardFunc,...
                                   InitEstFunc,...
                                   input_sigma,...
                                   params_admm,...
                                   orig_luma_im);
out_admm_im = MergeChannels(input_im,est_admm_im);

fprintf('Done.\n');


%% minimize the laplacian regularization functional via Steepest Descent

fprintf('Restoring using RED: Steepest-Descent method\n');

switch degradation_model
    case 'UniformBlur'
        params_sd = GetUniformDeblurSDParams(light_mode);
    case 'GaussianBlur'
        params_sd = GetGaussianDeblurSDParams(light_mode);
    case 'Downscale'
        params_sd = GetSuperResSDParams(light_mode);
    otherwise
        error('Degradation model is not defined');
end

[est_sd_im, psnr_sd] = RunSD(input_luma_im,...
                             ForwardFunc,...
                             BackwardFunc,...
                             InitEstFunc,...
                             input_sigma,...
                             params_sd,...
                             orig_luma_im);
% convert back to rgb if needed
out_sd_im = MergeChannels(input_im,est_sd_im);

fprintf('Done.\n');


%% display final results

fprintf('Image name %s \n', file_name);
fprintf('Input PSNR = %f \n', psnr_input);
fprintf('RED: Fixed-Point PSNR = %f \n', psnr_fp);
fprintf('RED: ADMM PSNR = %f \n', psnr_admm);
fprintf('RED: Steepest-Decent PSNR = %f \n', psnr_sd);


%% write images

if ~exist('./results/','dir')
    mkdir('./results/');
end

fprintf('Writing the images to ./results...');

imwrite(uint8(input_im),['./results/input_' file_name]);
imwrite(uint8(out_fp_im),['./results/est_fp_' file_name]);
imwrite(uint8(out_admm_im),['./results/est_admm_' file_name]);
imwrite(uint8(out_sd_im),['./results/est_sd_' file_name]);

fprintf(' Done.\n');

