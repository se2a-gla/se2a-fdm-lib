function blkStruct = slblocks
% This function specifies that the library should appear
% in the Library Browser
% and be cached in the browser repository

addpath(genpath(pwd));

Browser.Library = 'SE2A_library';
% 'mylib' is the name of the library

Browser.Name = 'SE2A library';
% 'My Library' is the library name that appears 
     % in the Library Browser

blkStruct.Browser = Browser; 