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

% This function creates a sparse matrix H that blurs (according to the
% input 'psf') and decimates (by factor 'scale') a high-resolution input
% image of size lr_im_sz*scale, resulting in a low-resolution image 
% of size lr_im_sz

% Note: This code is close to the function Set_blur_matrix, which 
% is provided in the NCSR package, given in:
% link: http://www4.comp.polyu.edu.hk/~cslzhang/NCSR.htm

function  H = CreateBlurAndDecimationOperator(scale, lr_im_sz, psf)

hr_im_sz = lr_im_sz*scale;

offset = (size(psf, 1)-1)/2;
cntr = ceil(size(psf, 1)/2);

tot_numel = size(psf, 1)^2 * lr_im_sz(1)*lr_im_sz(2);
row_locs = zeros(tot_numel, 1);
col_locs = zeros(tot_numel, 1);
values = zeros(tot_numel, 1);
cnt_numel = 1;

locs = 1 : 1 : (hr_im_sz(1)*hr_im_sz(2));
locs = reshape(locs, [hr_im_sz(1), hr_im_sz(2)]);

for lr_row = 1 : 1 : lr_im_sz(1)

    for lr_col = 1 : 1 : lr_im_sz(2)
        
        % find locations
        hr_row = (lr_row-1)*scale + 1;
        hr_col = (lr_col-1)*scale + 1;        
        
        row_idx = (lr_col-1)*lr_im_sz(1) + lr_row;
        rows_vec = max(hr_row - offset,1) : min(hr_row + offset,hr_im_sz(1));
        cols_vec = max(hr_col - offset,1) : min(hr_col + offset,hr_im_sz(2));
        col_idx = reshape(locs(rows_vec, cols_vec), [],1);
        cur_numel = size(col_idx,1);        
        row_locs(cnt_numel : cnt_numel+cur_numel-1) = row_idx;
        col_locs(cnt_numel : cnt_numel+cur_numel-1) = col_idx;
        
        % fill the values of the psf in the desired matrix
        rows_vec = cntr - (hr_row - max(hr_row - offset,1)) : cntr + min(hr_row + offset,hr_im_sz(1)) - hr_row;
        cols_vec = cntr - (hr_col - max(hr_col - offset,1)) : cntr + min(hr_col + offset,hr_im_sz(2)) - hr_col;
        crop_psf = reshape( psf(rows_vec,cols_vec), [], 1);
        values(cnt_numel : cnt_numel+cur_numel-1) = crop_psf/sum(crop_psf);
        
        cnt_numel = cnt_numel + cur_numel;

    end
    
end

row_locs = row_locs(1:cnt_numel-1);
col_locs = col_locs(1:cnt_numel-1);
values = values(1:cnt_numel-1);

H = sparse(...
    row_locs,...
    col_locs,...
    values,...
    lr_im_sz(1)*lr_im_sz(2),...
    hr_im_sz(1)*hr_im_sz(2));

return

