function tKey = bibkey
% bibkey Returns the key for the IEEE Xplore Metadata API
%
%   Usage:
%
%     bibkey
%
%   Description:
%
%     This purpose of this file is to contain the key for the
%     IEEE Xplore Metadata API, which can be obtained from:
%     https://developer.ieee.org/member/register
%
% see also bibget
% insert your key for the IEEE Xplore Metadata API below
tKey = '';
if isempty(tKey)
  fprintf('\n');
  fprintf('Each user of bibget requires a personal key for the\n');
  fprintf('IEEE Xplore Metadata API. You can get the key from\n');
  fprintf('\n');
  fprintf('  https://developer.ieee.org/member/register\n');
  fprintf('\n');
  fprintf('After registration, type ''open bibkey'' in the Matlab\n');
  fprintf('command window and store the key in the variable ''tKey''.\n');
  fprintf('\n');
  if nargout == 0
    clear tKey;
  else
    tKey = [];
  end
end
