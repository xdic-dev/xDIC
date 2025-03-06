function output = isText(myVar)
% isText checks if the input variable is a text type.
% The function returns true if the variable is of string, character array, or cell array of strings type; otherwise, it returns false.
%
% Input:
%   - myVar: The variable to be checked.
%
% Output:
%   - output: A logical value indicating if the variable is a text type.

    output = isstring(myVar) || ischar(myVar) || iscellstr(myVar);
end