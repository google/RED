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

% Objective:
%   Minimize E(x) = 1/(2sigma^2)||Hx-y||_2^2 + 0.5*lambda*x'*(x-denoise(x))
%   via steepest descend method.
%   Please refer to Section 4.2 in the paper for more details:
%   "Deploying the Denoising Engine for Solving Inverse Problems -- 
%   Gradient Descent Methods".
%
% Inputs:
%   y - the input image
%   ForwardFunc - the degradation operator H
%   BackwardFunc - the transpose of the degradation operator H
%   InitEstFunc - special initialization (e.g. the output of other method)
%   input_sigma - noise level
%   params.lambda - regularization parameter
%   params.outer_iters - number of total iterations
%   params.effective_sigma - the input noise level to the denoiser
%   orig_im - the original image, used for PSNR evaluation ONLY
%
% Outputs:
%   im_out - the reconstructed image
%   psnr_out - PSNR measured between x_est and orig_im

function [im_out, psnr_out] = RunSD(y, ForwardFunc, BackwardFunc,...
    InitEstFunc, input_sigma, params, orig_im)

% print info every PRINT_MOD steps
QUIET = 0;
PRINT_MOD = floor(params.outer_iters/10);
if ~QUIET
    fprintf('%7s\t%10s\t%12s\n', 'iter', 'PSNR', 'objective');
end

% parameters
lambda = params.lambda;
outer_iters = params.outer_iters;
effective_sigma = params.effective_sigma;

% compute step size
mu = 2/(1/(input_sigma^2) + lambda);

% initialization
x_est = InitEstFunc(y);

for k = 1:1:outer_iters
    
    % denoise
    f_x_est = Denoiser(x_est, effective_sigma);
    
    % update the solution
    grad1 = BackwardFunc(ForwardFunc(x_est) - y)/(input_sigma^2);
    grad2 = lambda*(x_est - f_x_est);
    x_est = x_est - mu*(grad1 + grad2);
    
    % project to [0,255]
    x_est = max( min(x_est, 255), 0);
    
    if ~QUIET && (mod(k,PRINT_MOD) == 0 || k == outer_iters)
        % evaluate the cost function
        fun_val = CostFunc(y, x_est, ForwardFunc, input_sigma,...
            lambda, effective_sigma);
        im_out = x_est(1:size(orig_im,1), 1:size(orig_im,2));
        psnr_out = ComputePSNR(orig_im, im_out);
        fprintf('%7i %12.5f %12.5f \n', k, psnr_out, fun_val);
    end
end

im_out = x_est(1:size(orig_im,1), 1:size(orig_im,2));
psnr_out = ComputePSNR(orig_im, im_out);

return

