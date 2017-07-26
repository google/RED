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

% Returns set of parameters
% set light_mode = true to run the code in a sub optimal but faster mode
% set light_mode = false to obtain the results reported in the RED paper

function params = GetSuperResSDParams(light_mode)

% regularization factor
params.lambda = 0.008;

% number of outer iterations
if light_mode
    params.outer_iters = 500;
else
    params.outer_iters = 1500;
end

% level of noise assumed in the regularization-denoiser
params.effective_sigma = 3;

return

