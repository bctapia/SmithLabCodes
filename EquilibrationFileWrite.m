% copies and changes the temperature for the equilibration script
clear; clc;

existquestion = input("Does a folder exist? Y/N [Y]: ","s");
if isempty(existquestion) || existquestion == 'Y'
    exist = 'Y';
elseif exist == 'N'
    disp('Not yet supported')
end

namequestion = input("Is the LAMMPS file named equil.in ? Y/N [Y]: ","s");
if isempty(namequestion) || namequestion == 'Y'
    filename = 'equil.in';
elseif exist == 'N'
    filename = input('file name: ', 's');
end

original_temperature = input('What is the existing folder temperature?: ');

% Prompt the user for desired temperatures
desired_temperatures = input('Enter desired temperatures separated by spaces: ', 's');
desired_temperatures = str2double(strsplit(desired_temperatures));

original_folder_path = string(original_temperature);
original_file_name = filename;
original_file_path = fullfile(original_folder_path, original_file_name);


% Loop through each desired temperature
for temp_idx = 1:numel(desired_temperatures)
    % Define the new folder path with temperature name
    new_temperature_folder_path = fullfile(num2str(desired_temperatures(temp_idx)));

    % Create the new folder
    mkdir(new_temperature_folder_path);

    % Copy all files from the original folder to the new folder
    copyfile(fullfile(original_folder_path, '*'), new_temperature_folder_path);

    % Define the new file name with temperature suffix
    new_file_name = sprintf('equil.in');

    % Define the new file path
    new_file_path = fullfile(new_temperature_folder_path, new_file_name);

    % Copy the original file to preserve it
    copyfile(original_file_path, new_file_path);

    % Open the original file for reading
    fileID = fopen(original_file_path, 'r');

    % Read the file line by line and store its contents
    file_contents = {};
    line = fgetl(fileID);
    while ischar(line)
        file_contents{end+1} = line;
        line = fgetl(fileID);
    end

    % Close the original file
    fclose(fileID);

    % Define the number to search for and replace
    new_temperature = desired_temperatures(temp_idx); % Use the current desired temperature

    % Loop through the file contents and replace the temperature
    for i = 1:numel(file_contents)
        % Find and replace the original temperature in each line
        file_contents{i} = strrep(file_contents{i}, num2str(original_temperature), num2str(new_temperature));
    end

    % Open the new file for writing
    fileID = fopen(new_file_path, 'w');

    % Write the modified contents to the new file
    for i = 1:numel(file_contents)
        fprintf(fileID, '%s\n', file_contents{i});
    end

    % Close the new file
    fclose(fileID);

    disp(['Replacement completed for temperature ', num2str(new_temperature)]);
end