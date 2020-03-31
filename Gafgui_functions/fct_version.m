
% --------------------------------------------------------------------
function err = fct_version()

%dumb = sprintf('Gafgui version 4.0\nCopyright (C) 2016-2020, by Hugo Bouchard.\nThis program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions. See code header for details.');
%msgbox(dumb,'Version');

%dumb = fileread('licenceGafgui.txt')
s = sprintf('Gafgui version 4.0');
s = sprintf('%s\nCopyright (C) 2007-2020 by Hugo Bouchard',s); 
s = sprintf('%s\nMy contact is h.bouchard@umontreal.ca',s);
s = sprintf('%s\n',s);  
s = sprintf('%s\nThis program is free software: you can redistribute it and/or modify',s); 
s = sprintf('%s\nit under the terms of the GNU General Public License as published by',s); 
s = sprintf('%s\nthe Free Software Foundation, either version 3 of the License, or',s); 
s = sprintf('%s\n(at your option) any later version.',s); 
s = sprintf('%s\n',s); 
s = sprintf('%s\nThis program is distributed in the hope that it will be useful,',s); 
s = sprintf('%s\nbut WITHOUT ANY WARRANTY; without even the implied warranty of',s); 
s = sprintf('%s\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the',s); 
s = sprintf('%s\nGNU General Public License for more details.',s); 
s = sprintf('%s\n',s);       
s = sprintf('%s\nYou should have received a copy of the GNU General Public License',s); 
s = sprintf('%s\nalong with this program.  If not, see <http://www.gnu.org/licenses/>.',s); 
s = sprintf('%s\n',s);  
s = sprintf('%s\nBy using the software, you agree to its conditions.',s); 

button = questdlg(s,'Version','I agree','I disagree','I agree');

if strcmp(button,'I agree')
    err = 0;
else
    err = 1;
end