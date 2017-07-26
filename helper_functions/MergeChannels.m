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

% Replaces the luminance channel of in_im with enhanced_im

function out_im = MergeChannels(in_im,enhanced_im)

if  size(in_im,3) == 3
    out_im = rgb2ycbcr(uint8(in_im));   
    out_im(:,:,1) = enhanced_im;
    out_im = ycbcr2rgb(uint8(out_im));   
else
    out_im = enhanced_im;
end

