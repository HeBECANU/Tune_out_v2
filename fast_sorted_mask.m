function idxs=fast_sorted_mask(data,min_val,max_val)
%fast_sorted_mask - a fast replacement for masking of ordered vectors based on binary search
% DATA VECTOR MUST BE SORTED! THERE IS NO CHECK, IT WILL JUST RETURN THE WRONG ANSWER!!!
% gives a speedup (~100) for taking small slices of large (>1e6 elments) sorted vectors
% if the data is not already sorted then you can still get a net speedup
% by sorting the data then use this code upwards of 30 times 
% see test_fast_sorted_mask for a number of speed comparisons

% Syntax:  min_max_val_idx=fast_sorted_mask(data,min_val,max_val)
% Example: 
%       %show speedup, basic masking 19ms this code 0.5ms
%       data=sort(rand(1e7,1));
%       min_val=0.9;
%       max_val=0.91;
%       tic; mask=data<max_val & data>min_val;
%       subdata1=sort(data(mask)); toc
%       tic;mask_idx=fast_sorted_mask(data,min_val,max_val);
%       subdata2=data(mask_idx(1):mask_idx(2)); toc
%       isequal(subdata1,subdata2)

% Other m-files required: none
% Also See:test_fast_sorted_mask
% Subfunctions: binary_search_first_elm
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
% 	-fix whatever is causeing execution time to flatten out in log-log past ~7e6.  see test_fast_sorted_mask
% 	-should compile it for more speed!

% Author: Bryce Henson
% email: Bryce.Henson[a circle]live.com  %YOU MUST INCLUDE '[matlab][fast_sorted_mask]' in the subject line OR I WILL NOT REPLY
% Last revision:2018-09-30

%------------- BEGIN CODE --------------


elms=numel(data);

if data(1)>min_val
    min_idx=1;
    if data(1)>max_val
        max_idx=min_idx-1;%no output
    end
else
    min_idx=binary_search_first_elm(data,min_val,1,elms);
    if data(min_idx)<min_val
        min_idx=min_idx+1;
    end
end
   
if data(end)<max_val
    max_idx=elms;
    if data(end)<min_val
        min_idx=max_idx+1;%no output
    end
else
    max_idx=binary_search_first_elm(data,max_val,min_idx,elms);
    if data(max_idx)>max_val
        max_idx=max_idx-1;
    end
end

% if max_idx<min_idx
%     warning('fast_sorted_mask:no data points in range')
% end

idxs=[min_idx,max_idx];

end


%modified from mathworks submission by Benjamin Bernard 
%from https://au.mathworks.com/matlabcentral/fileexchange/37915-binary-search-for-closest-value-in-an-array
function idx = binary_search_first_elm(vec, val,min_idx,start_idx)
% Returns index of vec that is closest to val, searching between min_idx start_idx . 
%If several entries
% are equally close, return the first. Works fine up to machine error (e.g.
% [v, i] = closest_value([4.8, 5], 4.9) will return [5, 2], since in float
% representation 4.9 is strictly closer to 5 than 4.8).
% ===============
% Parameter list:
% ===============
% arr : increasingly ordered array
% val : scalar in R

top = start_idx(1);
btm = min_idx(1);

% Binary search for index
while top - btm > 1
    med = floor((top + btm)/2);
    % Replace >= here with > to obtain the last index instead of the first.
    if vec(med) >= val 
        top = med;
    else
        btm = med;
    end
end

idx=btm;
end

% Copyright (c) 2012, Benjamin Bernard
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
