function write_file(mat_string, file_path)
% Writes a cell array or matrix to a text file.
%
% Input:
%   - mat_string: Cell array or matrix to be written to the file.
%   - file_path: File path for writing the content.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = size(mat_string, 1);
    fid = fopen(file_path, 'wt');
    
    for i = 1:n
        if iscell(mat_string)
            current_line = mat_string{i};
        else
            current_line = mat_string(i, :);
        end
        fprintf(fid, '%c', current_line);
        fprintf(fid, '\n');
    end
    
    fclose(fid);
end
