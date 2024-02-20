function cleaned_str = only_letters(input_str, varargin)
% Replaces non-letter characters in input_str with spaces.
%
% Input:
%   - input_str: Input string or cell array of strings.
%   - varargin: Optional argument to include numbers (if set to 1).
%
% Output:
%   - cleaned_str: Cleaned string or cell array of strings.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    V = [];
    if ~isempty(varargin)
        V = varargin{1};
    end
    
    if iscell(input_str)
        cleaned_str = trim_all(cellfun(@(x) only_letters_str(x, V), ...
            input_str, 'UniformOutput', false));
    else
        cleaned_str = only_letters_str(input_str, V);
    end
end

function cleaned_str = only_letters_str(input_str, varargin)
% Replaces non-letter characters in input_str with spaces.
    put_num_off = [];
    
    if ~isempty(varargin)
        if varargin{1} == 1
            put_num_off = 48:57; % ASCII values for numbers
        end
    end
    
    % Characters to be replaced with space
    not_letters = [247 215 216 123:191 91:96 1:47 58:64 put_num_off];
    
    in_list = ismember(input_str, not_letters);
    cleaned_str = input_str;
    cleaned_str(in_list) = ' ';
end
