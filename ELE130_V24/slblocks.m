function blkStruct = slblocks
% This function specifies that the library 'mylib'
% should appear in the Library Browser with the
% name 'Eget bibliotek'

Browser.Library = 'mylib';
% 'mylib' is the name of the library

Browser.Name = 'Eget bibliotek';
% 'Eget bibliotek' is the library name that appears
% in the Library Browser

blkStruct.Browser = Browser;