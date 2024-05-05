function process_lmps_file(input_file, output_file)
    % Open the original file
    fileID = fopen(input_file, 'r');

    % Read all lines from the file
    all_lines = textscan(fileID, '%s', 'Delimiter', '\n');
    all_lines = all_lines{1};

    % Close the original file
    fclose(fileID);

    % Find the index of the line containing "Atoms # full"
    atoms_index = find(strcmp(all_lines, 'Atoms # full'));

    % Initialize the index for data extraction
    data_start_index = atoms_index + 1;

    % Find the index of the second blank line after the data following "Atoms # full"
    data_end_index = data_start_index;
    blank_line_count = 0;
    while blank_line_count < 2 && data_end_index <= numel(all_lines)
        if isempty(all_lines{data_end_index})
            blank_line_count = blank_line_count + 1;
        end
        data_end_index = data_end_index + 1;
    end
    data_end_index = data_end_index - 1;

    % Extract the data lines between "Atoms # full" and "Velocities"
    data_lines = all_lines(data_start_index:data_end_index);

    % Parse the data into a matrix
    data_matrix = [];
    for i = 1:numel(data_lines)
        % Split each line by spaces and convert to double
        row = str2double(strsplit(data_lines{i}, {' ', '\t'}));
        % Remove any NaN entries
        row = row(~isnan(row));
        data_matrix = [data_matrix; row];
    end

    % Sort the data matrix based on the first column
    sorted_data_matrix = sortrows(data_matrix, 1);

    % Open a new file to write the updated data
    new_fileID = fopen(output_file, 'w');

    % Write lines before "Atoms # full" unchanged
    for i = 1:atoms_index
        fprintf(new_fileID, '%s\n', all_lines{i});
    end

    % Write a newline after "Atoms # full"
    fprintf(new_fileID, '\n');

    % Write the sorted data lines with formatted values
    for i = 1:size(sorted_data_matrix, 1)
        line = sorted_data_matrix(i, :);
        formatted_line = '';
        for j = 1:numel(line)
            value = line(j);
            if mod(value, 1) == 0 % Check if the number is an integer
                formatted_line = [formatted_line, sprintf('%d', value)];
            else
                formatted_line = [formatted_line, sprintf('%.6f', value)];
            end
            if j < numel(line)
                formatted_line = [formatted_line, ' ']; % Add space delimiter
            end
        end
        fprintf(new_fileID, "%s\n", formatted_line);
    end

    % Write a blank line before "Velocities"
    fprintf(new_fileID, '\n');

    % Write "Velocities"
    fprintf(new_fileID, 'Velocities\n');

    % Write a blank line after "Velocities"
    fprintf(new_fileID, '\n');

    % Find the index of "Velocities" and "Bonds"
    velocities_index = find(strcmp(all_lines, 'Velocities'));
    bonds_index = find(strcmp(all_lines, 'Bonds'));

    % Extract the data lines between "Velocities" and "Bonds"
    velocities_lines = all_lines(velocities_index+1:bonds_index-1);

    % Parse the velocities data into a matrix
    velocities_matrix = [];
    for i = 1:numel(velocities_lines)
        % Split each line by spaces and convert to double
        row = str2double(strsplit(velocities_lines{i}, {' ', '\t'}));
        % Remove any NaN entries
        row = row(~isnan(row));
        velocities_matrix = [velocities_matrix; row];
    end

    % Sort the velocities matrix based on the first column
    sorted_velocities_matrix = sortrows(velocities_matrix, 1);

    % Write the sorted velocities data to the file
    for i = 1:size(sorted_velocities_matrix, 1)
        line = sorted_velocities_matrix(i, :);
        formatted_line = sprintf('%d ', line);
        fprintf(new_fileID, "%s\n", formatted_line);
    end

    % Write a blank line before "Bonds"
    fprintf(new_fileID, '\n');

    % Write "Bonds"
    fprintf(new_fileID, 'Bonds\n');

    % Write the lines after "Bonds" unchanged until the end of the file
    for i = bonds_index+1:numel(all_lines)
        fprintf(new_fileID, '%s\n', all_lines{i});
    end

    % Close the new file
    fclose(new_fileID);
    disp(['File updated: ' output_file]); % Print notification
end


function process_files_in_folder(folder_path)
    % List all files and folders in the folder
    contents = dir(folder_path);
    
    % Iterate over each item
    for i = 1:length(contents)
        item = contents(i);
        if item.isdir && ~strcmp(item.name, '.') && ~strcmp(item.name, '..')
            % If the item is a folder (and not '.' or '..'), recursively call the function
            subfolder_path = fullfile(folder_path, item.name);
            process_files_in_folder(subfolder_path);
        elseif ~item.isdir && startsWith(item.name, 'CANAL')
            % If the item is a file and its name starts with 'CANAL', process it
            file_path = fullfile(folder_path, item.name);
            
            % Create a new file in the same folder
            [~, base_name, ext] = fileparts(item.name);
            new_file_name = [base_name '_new' ext];
            new_file_path = fullfile(folder_path, new_file_name);
            fid = fopen(new_file_path, 'w');
            fclose(fid);
            
            % Call your function on this file with the new file path as the second argument
            process_lmps_file(file_path, new_file_path);
        end
    end
end

% Base folder path
base_folder_path = 'Folder'; % Replace this with the path to your base folder

% Start processing files in the base folder
process_files_in_folder(base_folder_path);

exit
