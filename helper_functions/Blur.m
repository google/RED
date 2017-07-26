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

% Blurs the image x with an input psf. This process is similar to the 
% degradation flow of the NCSR paper
function out_im = Blur(in_im, psf)

cntr = (size(psf)-1)/2;

% treat the borders
pad_im = [in_im(:,cntr:-1:1,:), in_im, in_im(:,end:-1:end-cntr+1,:)];
pad_im = [pad_im(cntr:-1:1,:,:); pad_im; pad_im(end:-1:end-cntr+1,:,:)];

if size(in_im,3) == 3
    out_im(:,:,1) = conv2(pad_im(:,:,1), psf, 'valid');
    out_im(:,:,2) = conv2(pad_im(:,:,2), psf, 'valid');
    out_im(:,:,3) = conv2(pad_im(:,:,3), psf, 'valid');
else
    out_im = conv2(pad_im, psf, 'valid');
end
