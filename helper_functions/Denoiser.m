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

% Calls a denoiser and return the cleaned image f_est.
% In this example the chosen denoising method is TNRD.
% Note: The devision and multiplication by the constant 5 is done to handle
% a different noise-level than a fixed one (which is equal to 5)

function f_est = Denoiser(x_est,effective_sigma)

f_est = ReactionDiffusion(5/effective_sigma*x_est);
f_est = f_est*effective_sigma/5;

return

