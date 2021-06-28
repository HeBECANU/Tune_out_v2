function bibget(tBibKey, varargin)
% bibget The easiest way to get BibTeX entries from IEEE Xplore
%
%   Usage:
%
%     bibget NameYYx
%
%     bibget NameYYx Keyword1 Keyword2 ...
%
%     bibget NameYYx Option Keyword1 Keyword2 ...
%
%     Option        Description
%     ---------------------------------------------------------------------
%     -nosplit      Do not split camelcase names for search.
%                   Search for "McDonald", instead of "Mc Donald".
%     -overwrite    Send a new query to IEEE Xplore and
%                   overwrite the current entry in the BibTeX file.
%     -NUMBER       Directly identify an article by its article number.
%                   NUMBER is a placeholder for an article number.
%
%   Description:
%
%     The function bibget takes one or more string input arguments,
%     specifying a bibliographic search in the IEEE Xplore database.
%     After querying the IEEE Xplore database, it lists the results
%     in the command window and asks the user to select one result.
%     The infos of the selected result are converted to a BibTeX entry,
%     which is displayed and stored in a BibTeX file.
%
%     By default, bibget uses the BibTeX file in the current working
%     directory of Matlab. If there is no BibTeX file, it creates one
%     named 'bib.bib'. If there are multiple BibTeX files in the current
%     working directory, bibget asks the user to select one of them.
%
%     Additional features of bibget are:
%      - Conversion of special characters in the title to Unicode
%      - Protection of acronyms in the title by curly braces
%      - Removal of numeric infos from the publication title
%     
%   Syntax:
%
%     The first input argument must have the form NameYYx, where Name is
%     the last name of the first author, YY are the last two digits of
%     the publication year and x is an optional lowercase letter to allow
%     the storing of more than one publication per author and year.
%
%     If the last name is very common, it is likely that there are too
%     many results to be displayed. For this reason it is possible to
%     specify additional search keywords as optional input arguments.
%
%     If the last name consists of more than one part and the parts
%     are separated by spaces, dashes (-) or primes (') it should be
%     written in camelcase i.e. FirstSecondThird.
%
%     If the last name includes camelcase that should not be separated
%     like in McDonald, the -nosplit option should be used.
%
%     If the last name includes typographic accents or other non-ASCII
%     characters, a transliteration to ASCII should be used.
%
% see also bibload, bibsave
% written by      Harald Enzinger
% version         2.4
% released at     www.mathworks.com/matlabcentral/fileexchange/53412
% release date    01.03.2019
% key for IEEE Xplore Metadata API
tApiKey = bibkey;
% stop if there is no key
if isempty(tApiKey)
  return;
end
% check first input argument
if nargin < 1
  fprintf('\n');
  fprintf('There must be at least one input argument.\n');
  fprintf('\n');
  return;
elseif ~ischar(tBibKey) || (size(tBibKey,2) ~= numel(tBibKey))
  fprintf('\n');
  fprintf('The first input argument must be a string.\n');
  fprintf('\n');
  return;
end
% initialize variables for command line options
bNoSplit    = false;
bOverwrite  = false;
bArNumber   = false;
tArNumber   = [];
% check additional input arguments
if nargin > 1
  cKeywords = varargin;
  nKeywords = length(cKeywords);
  vKeywords = false(1,nKeywords);
  for iKeyword = 1 : nKeywords
    % check data type
    if ischar(cKeywords{iKeyword})
      tKeyword = cKeywords{iKeyword};
    else
      fprintf('\n');
      fprintf('All optional input arguments must be strings.\n');
      fprintf('\n');
      return;
    end
    % check semantic type
    if tKeyword(1) == '-'
      cTok = regexp(tKeyword,'-(\d+)$','tokens');
      if ~isempty(cTok)
        bArNumber = true;
        tArNumber = cTok{1}{1};
      else
        switch(tKeyword)
          case '-nosplit'
            bNoSplit = true;
          case '-overwrite'
            bOverwrite = true;
          otherwise
            fprintf('\n');
            fprintf('Unknown option string %s.\n', tKeyword);
            fprintf('\n');
            return;
        end
      end
    else
      vKeywords(iKeyword) = true;
    end
  end
  % extract keywords
  cKeywords = cKeywords(vKeywords);  
else
  cKeywords = {};
end
% parse bibkey
rBibKey = parseBibKey(tBibKey, bNoSplit);
% stop here if rBibKey is empty (i.e. wrong format of tBibKey)
if isempty(rBibKey)
  return;
end
% get bibtex data structure (this loads or creates a new file)
rBib = bibload;
% stop here if rBib is empty (i.e. aborted by the user)
if isempty(rBib)
  return;
end
% check if bibtex key is already present in the data structure
[bMember, iMember] = ismember(tBibKey, rBib.cKey);
if bMember && ~bOverwrite
  fprintf('\n');
  disp(rBib.cBib{iMember});
  fprintf('\n');
  return;
end
% get article by number or search for it
if bArNumber
  hDoc = getArticle(tApiKey, tArNumber, bOverwrite);
else
  hDoc = searchArticle(tApiKey, rBibKey, cKeywords, bOverwrite);
end
% stop here if no article was found
if isempty(hDoc)
  return;
end
% get author string
hAuthors = hDoc.getElementsByTagName('full_name');
nAuthors = hAuthors.getLength;
cAuthors = cell(2,nAuthors);
for iAuthor = 1 : nAuthors
  cAuthors{1,iAuthor} = char(hAuthors.item(iAuthor-1).item(0).getData);
  if iAuthor == nAuthors
    cAuthors{2,iAuthor} = '';
  else
    cAuthors{2,iAuthor} = ' and ';
  end
end
tAuthors = [cAuthors{:}];
% get other info strings
tPubType    = getData(hDoc, 'content_type');
tTitle      = getData(hDoc, 'title');
tPubTitle   = getData(hDoc, 'publication_title');
tVolume     = getData(hDoc, 'volume');              % journal papers only
tNumber     = getData(hDoc, 'issue');               % journal papers only
tLocation   = getData(hDoc, 'conference_location'); % conference papers only
tPageStart  = getData(hDoc, 'start_page');
tPageEnd    = getData(hDoc, 'end_page');
tArNumber   = getData(hDoc, 'article_number');
if strcmp(tPubType, 'Conferences')
  tDate = getData(hDoc, 'conference_dates');
else
  tDate = getData(hDoc, 'publication_date');
end
% determine year and month
if isempty(tDate)
  % use year from bibkey
  tYear   = rBibKey.tYearQuery;
  tMonth  = '';
else  
  % extract year from date string
  cYear = regexp(tDate, '(19|20)[0-9]{2}', 'match');
  if isempty(cYear)
    tYear = '';
  else
    tYear = cYear{1};
  end
  % extract month from date string
  tMonths = 'Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec';
  cMonth  = regexp(tDate, tMonths, 'match');
  if isempty(cMonth)
    tMonth = '';
  else
    tMonth = cMonth{1};
    switch(tMonth)
      case 'Jan'; tMonth = 'January';
      case 'Feb'; tMonth = 'February';
      case 'Mar'; tMonth = 'March';
      case 'Apr'; tMonth = 'April';
      case 'Jun'; tMonth = 'June';
      case 'Jul'; tMonth = 'July';
      case 'Aug'; tMonth = 'August';
      case 'Sep'; tMonth = 'September';
      case 'Oct'; tMonth = 'October';
      case 'Nov'; tMonth = 'November';
      case 'Dec'; tMonth = 'December';
    end
  end    
end
% process title string
if ~isempty(tTitle)
  % -----------------------------------------------------------------------
  % remove html tags
  % -----------------------------------------------------------------------
  % remove spaces between html tags
  tTitle = regexprep(tTitle, '> <', '><');
  % replace img html tags by their alt attributes in math mode
  tTitle = regexprep(tTitle, '<img[^>]+alt="([^"]+)">', '$$1$');  
  % remove other html tags
  tTitle = regexprep(tTitle, '<[^>]+>', '');
  % -----------------------------------------------------------------------
  % unicode conversion
  % -----------------------------------------------------------------------
  tTitle = unicode(tTitle);
  % -----------------------------------------------------------------------
  % add backslash for protection of latex special characters
  % -----------------------------------------------------------------------
  tTitle = regexprep(tTitle, '([#$%&~_^\{}])', '\\$1');  
  % -----------------------------------------------------------------------
  % add curly braces for protection of meaningful uppercase
  % -----------------------------------------------------------------------
  % protect acronyms
  % [A-Z]-\S*          Q-value
  % \S*-[A-Z]          class-A
  % [^a-z\s]+[A-Z]     ADC I/Q 64-QAM 5G
  % [\w.]+[A-Z]+\w*    12.3dBm
  tTitle  = regexprep(tTitle, '\<([A-Z]-\S*|\S*-[A-Z]|[^a-z\s]+[A-Z]|[\w.]+[A-Z]+\w*)\>', '{$1}');
  % protect uppercase words following a colon
  tTitle  = regexprep(tTitle, ': ([A-Z][a-z]*)', ': {$1}');
  % create the character em dash
  tEmDash = char(8212);
  % enforce space around em dash
  tExpr   = ['(\S)',tEmDash,'(\S)'];
  tRepl   = ['$1 ',tEmDash,' $2'];
  tTitle  = regexprep(tTitle, tExpr, tRepl);
  % protect uppercase words following a space-separated hyphen or em dash
  tExpr   = ['( - | ',tEmDash,' )([A-Z][a-z]*)'];
  tTitle  = regexprep(tTitle, tExpr, '$1{$2}');
end
% process pubtitle string
if ~isempty(tPubTitle)  
  if strcmp(tPubType, 'Conferences')
    % remove year from beginning of conference title
    tPubTitle = regexprep(tPubTitle, ['^',tYear,' '], '');
  end
  % add backslash for protection of latex special characters
  tPubTitle = regexprep(tPubTitle, '([#$%&~_^\{}])', '\\$1');  
end
% combine page numbers to page string
if ~isempty(tPageStart) && ~isempty(tPageEnd)
  tPages = [tPageStart, '-', tPageEnd];
else
  tPages = [];
end
% initialize strings
tBib = '';
tSep = sprintf(',\n');
% construct bibtex entry
switch(tPubType)
  
  case {'Journals', 'Journals & Magazines', 'Early Access'}
    
    tBib = [tBib, '@ARTICLE{',          tBibKey];
    tBib = [tBib, tSep, 'author={',     tAuthors,     '}'];
    tBib = [tBib, tSep, 'title={',      tTitle,       '}'];
    tBib = [tBib, tSep, 'journal={',    tPubTitle,    '}'];
    tBib = [tBib, tSep, 'year={',       tYear,        '}'];
    tBib = [tBib, tSep, 'month={',      tMonth,       '}'];
    tBib = [tBib, tSep, 'volume={',     tVolume,      '}'];    
    tBib = [tBib, tSep, 'number={',     tNumber,      '}'];
    tBib = [tBib, tSep, 'pages={',      tPages,       '}'];    
    tBib = [tBib, tSep, 'arnumber={',   tArNumber,    '}'];
    tBib = [tBib, '}'];
    
  case 'Conferences'
    
    tBib = [tBib, '@INPROCEEDINGS{',    tBibKey];
    tBib = [tBib, tSep, 'author={',     tAuthors,     '}'];
    tBib = [tBib, tSep, 'title={',      tTitle,       '}'];
    tBib = [tBib, tSep, 'booktitle={',  tPubTitle,    '}'];
    tBib = [tBib, tSep, 'year={',       tYear,        '}'];      
    tBib = [tBib, tSep, 'month={',      tMonth,       '}'];      
    tBib = [tBib, tSep, 'pages={',      tPages,       '}'];      
    tBib = [tBib, tSep, 'location={',   tLocation,    '}'];      
    tBib = [tBib, tSep, 'arnumber={',   tArNumber,    '}'];      
    tBib = [tBib, '}'];
    
  otherwise
    
    fprintf('\n');
    fprintf('The publication type ''%s'' is unknown.\n', tPubType);
    fprintf('No bibtex entry was saved.\n');
    fprintf('\n');
    return;
    
end
% update or add bibtex entry
if bOverwrite && bMember
  rBib.cKey{iMember} = tBibKey;
  rBib.cBib{iMember} = tBib;
else
  rBib.cKey = [rBib.cKey; {tBibKey}];
  rBib.cBib = [rBib.cBib; {tBib}];
end
% sort bibtex entries according to bibtex key
[rBib.cKey, vIdx] = sort(rBib.cKey);
rBib.cBib = rBib.cBib(vIdx);
% update date number
rBib.tDate = datestr(now);
% assign bibtex data structure in base workspace
assignin('base', 'bib', rBib);
% save bibtex data structure to file
bibsave;
% display bibtex entry in command window
fprintf('\n');
disp(tBib);
fprintf('\n');
% --- Subfunction ---------------------------------------------------------
function rBibKey = parseBibKey(tBibKey, bNoSplit)
% parse first input argument
cBibKey = regexp(tBibKey, '^([A-Za-z]{2,})([0-9]{2})(?:[a-z]?)$', 'tokens');
if isempty(cBibKey)
  fprintf('\n');
  fprintf('The first input argument must have the form NameYYx,\n');
  fprintf('where Name is the last name of the first author,\n');
  fprintf('YY are the last two digits of the publication year,\n');
  fprintf('and x is an optional lowercase letter.\n');
  fprintf('\n');
  rBibKey = [];
  return;
else
  tName = cBibKey{1}{1};
  tYear = cBibKey{1}{2};
end
% determine camelcase flag
bCamelCase = ~isempty(regexp(tName, '.+[A-Z]', 'once'));
% determine author name for search
if bNoSplit
  % use name from bibtex-key without modification
  tNameQuery = tName;
else
  % insert url-spaces (%20) into camelcase
  tNameQuery = regexprep(tName, '(.)([A-Z])', '$1%20$2');
end
% publication year for search
tDate = date;
nYear12Now = str2double(tDate(8:9));
nYear34Now = str2double(tDate(10:11));
nYear34Key = str2double(tYear);
if nYear34Key <= nYear34Now
  nYear = nYear12Now * 100 + nYear34Key;
else
  nYear = (nYear12Now-1) * 100 + nYear34Key;
end
tYearQuery = num2str(nYear);
% store infos
rBibKey.tBibKey     = tBibKey;
rBibKey.tName       = tName;
rBibKey.tYear       = tYear;
rBibKey.bCamelCase  = bCamelCase;
rBibKey.bNoSplit    = bNoSplit;
rBibKey.tNameQuery  = tNameQuery;
rBibKey.tYearQuery  = tYearQuery;
% --- Subfunction ---------------------------------------------------------
function hDoc = getArticle(tApiKey, tArNumber, bOverwrite)
% construct url
tUrl = [  'http://ieeexploreapi.ieee.org/api/v1/search/articles', ...
          '?apikey=', tApiKey, ...
          '&format=xml', ...
          '&article_number=', tArNumber ];
% query IEEE Xplore
[hRes, nTot] = queryIeeeXplore(tUrl, bOverwrite);
% check number of results
if nTot == -1
  hDoc = [];
  return;
elseif nTot == 0
  fprintf('\n');
  fprintf('There is no article with article number %s.\n', tArNumber);
  fprintf('\n');
  hDoc = [];
  return;
end
% get article
hDoc = hRes.getElementsByTagName('article').item(0);
% --- Subfunction ---------------------------------------------------------
function hDoc = searchArticle(tApiKey, rBibKey, cKeywords, bOverwrite)
% construct url
tUrl = [  'http://ieeexploreapi.ieee.org/api/v1/search/articles', ...
          '?apikey=', tApiKey, ...
          '&format=xml', ...
          '&author=', rBibKey.tNameQuery, ...
          '&publication_year=', rBibKey.tYearQuery ];
% add keywords to url
if ~isempty(cKeywords)
  cKeywords   = [ cKeywords; ...
                  repmat({'%20'},[1,length(cKeywords)-1]), {''} ];
  tUrl        = [ tUrl, '&meta_data=', cKeywords{:} ];
end
% query IEEE Xplore
[hRes, nTot] = queryIeeeXplore(tUrl, bOverwrite);
% check number of results
if nTot == -1
  % message was already displayed
  hDoc = [];
  return;
elseif nTot == 0
  % display message
  if rBibKey.bCamelCase && ~rBibKey.bNoSplit
    fprintf('\n');
    fprintf('The search returned no results. (the -nosplit option might help)\n');
    fprintf('\n');
  else
    fprintf('\n');
    fprintf('The search returned no results.\n');
    fprintf('\n');
  end
  hDoc = [];
  return;
end
% abort if number of papers is too high
if nTot > 25
  fprintf('\n');
  fprintf('There are %d results. Please use additional keywords.\n', nTot);
  fprintf('\n');
  hDoc = [];
  return;
end
% use last word in name for matching
cNameMatch = regexp(rBibKey.tNameQuery, '[A-Z][a-z]+$', 'match');
tNameMatch = cNameMatch{1};
nNameMatch = length(tNameMatch);
% remove all papers where first-author does not match
hDocs = hRes.getElementsByTagName('article');
hAuth = hRes.getElementsByTagName('authors');
iItem = 0;
nDocs = hDocs.getLength;
for iDoc = 1 : nDocs
  % get author string
  nAuth = hAuth.item(iItem).getLength;
  if nAuth > 0
    hFullName = hAuth.item(0).getElementsByTagName('full_name');
    tAuthor   = char(hFullName.item(0).item(0).getData);
  else
    iItem = iItem + 1;
    continue;
  end
  % compare author string
  nChars  = length(tAuthor);
  nMatch  = min(nChars,nNameMatch);
  vAuthor = double(tAuthor(nChars:-1:nChars-nMatch+1));
  vMatch  = double(tNameMatch(nNameMatch:-1:nNameMatch-nMatch+1));  
  for iChar = 1 : nMatch
    if vAuthor(iChar) < 128 && vAuthor(iChar) ~= vMatch(iChar)
      % remove article
      hDoc = hDocs.item(iItem);
      hRes.item(0).removeChild(hDoc);
      iItem = iItem - 1;  
      break;
    end
  end
  iItem = iItem + 1;
end
% remaining first-author papers
hTitles = hRes.getElementsByTagName('title');
nTitles = hTitles.getLength;
% abort if no first-author papers were found
if nTitles == 0
  fprintf('\n');
  fprintf('No first-author papers were found.\n');
  fprintf('\n');
  hDoc = [];
  return;
end
% display titles  
fprintf('\n');
for iTitle = 1 : nTitles
  % get title
  if hTitles.item(iTitle-1).getLength > 0
    tTitle = char(hTitles.item(iTitle-1).item(0).getData);
  else
    tTitle = 'title not available';
  end
  % remove all html tags
  tTitle = regexprep(tTitle, '<[^>]+>', '');
  % remove multiple white-spaces
  tTitle = regexprep(tTitle, '\s+', ' ');
  % remove white-spaces from begin or end of line
  tTitle = regexprep(tTitle, '^\s+|\s+$', '');
  % replace unicode references by unicode characters
  tTitle = unicode(tTitle);
  % display line
  fprintf(' %2d - %s\n', iTitle, tTitle);
end
fprintf('\n');
% request user input
while true
  tSel = input('Please select a number: ', 's');
  if isempty(tSel)
    fprintf('\baborted by the user\n');
    fprintf('\n');
    hDoc = [];
    return;
  else
    iSel = str2double(tSel);
    if isnan(iSel) || mod(iSel,1) || iSel < 1 || iSel > nTitles
      fprintf('\b invalid input!\n');
    else
      break;
    end
  end
end
% get article
hDoc = hRes.getElementsByTagName('article').item(iSel-1);
% --- Subfunction ---------------------------------------------------------
function [hRes, nTot] = queryIeeeXplore(tUrl, bOverwrite)
% path to temporary files
tFileUrl = [tempdir, 'getbib.url'];
tFileXml = [tempdir, 'getbib.xml'];
% new query flag
if bOverwrite || ~exist(tFileUrl, 'file') || ~exist(tFileXml, 'file')
  bNewQuery = true;
else
  % get url from previous call
  hFile = fopen(tFileUrl, 'r');
  tUrlPrev = fgetl(hFile);
  fclose(hFile);
  % compare new and previous url
  bNewQuery = ~strcmp(tUrl, tUrlPrev);
end
% perform new query
if bNewQuery
  try
    % send query and store result in xml file  
    websave(tFileXml, tUrl);
  catch oException
    % display message
    switch(oException.identifier)
      case 'MATLAB:webservices:HTTP403StatusCodeError'
        fprintf('\n');
        fprintf('The server denied your request.\n');
        fprintf('\n');
        fprintf('Either your API key is not yet activated, or\n');
        fprintf('you reached the maximum of 200 requests per day.\n');
        fprintf('\n');
      case 'MATLAB:webservices:HTTP504StatusCodeError'
        fprintf('\n');
        fprintf('A timeout occured. Please try again.\n');
        fprintf('\n');        
      otherwise
        oException.throw;
    end
    hRes = [];
    nTot = -1;
    % store empty string in url file
    hFile = fopen(tFileUrl, 'w');
    fprintf(hFile, '');
    fclose(hFile);    
    return;
  end
  % store url in url file
  hFile = fopen(tFileUrl, 'w');
  fprintf(hFile, '%s', tUrl);
  fclose(hFile);  
end
% parse xml file
try
  hRes = xmlread(tFileXml);
catch
  fprintf('\n');
  fprintf('Unable to read the file returned by IEEE Xplore.\n');
  fprintf('Query URL: %s\n', tUrl);
  fprintf('\n');
  hRes = [];
  nTot = -1;
  return;       
end
% total number of papers found
hTot = hRes.getElementsByTagName('totalfound');
if hTot.getLength
  nTot = str2double(char(hTot.item(0).item(0).getData));
else
  nTot = 0;
end
% --- Subfunction ---------------------------------------------------------
function tData = getData(hDoc, tField)
% get data field from xml structure
hField = hDoc.getElementsByTagName(tField);
if hField.getLength > 0 && hField.item(0).getLength > 0
  tData = char(hField.item(0).item(0).getData);
else
  tData = [];
end
