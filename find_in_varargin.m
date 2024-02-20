function value = find_in_varargin(varargin_input, target_name, default_value)
% Find in Varargin
%
% This function returns the value assigned to the target_name if it occurs in
% varargin_input (external varargin).
%
% Parameters:
%   varargin_input: External varargin cell array
%   target_name: Name of the parameter to find in varargin
%   default_value: Default value to return if target_name is not found
%
% Output:
%   value: Value assigned to target_name if found, otherwise default_value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the default value
value = default_value;

% Check if varargin_input is not empty
if ~isempty(varargin_input)
    % Loop through varargin pairs and find the value assigned to target_name
    for ii = 1:length(varargin_input)
        current_element = varargin_input{ii};
        
        % Check if the current element is a character and matches the target_name
        if ischar(current_element) && strcmpi(current_element, target_name)
            % Update the value if a match is found
            value = varargin_input{ii + 1};
        end
    end
end

end
