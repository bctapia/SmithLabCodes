clear; clc; close all

s = 4;
f = 20003;
gas = 'CH4';
temp = '308';
Trial = '3';

path = strcat('\\wsl.localhost\Ubuntu\home\btapia\MCMD\', temp, '\', gas, '\Trial', Trial, '_', gas);

% 0.1 %
iters = 0;
result01 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result01(:, i) = MCMD_import(strcat(path, '\0.1\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 0.5 %
iters = 0;
result05 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result05(:, i) = MCMD_import(strcat(path, '\0.5\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 1 %
iters = 10;
result1 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result1(:, i) = MCMD_import(strcat(path, '\1\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 5 %
iters = 10;
result5 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result5(:, i) = MCMD_import(strcat(path, '\5\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 10 %
iters = 0;
result10 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result10(:, i) = MCMD_import(strcat(path, '\10\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 20 %
iters = 10;
result20 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result20(:, i) = MCMD_import(strcat(path, '\20\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 30 %
iters = 10;
result30 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result30(:, i) = MCMD_import(strcat(path, '\30\results\', num2str(i), '.gcmc.prp'), s, f);
end

% 40 %
iters = 10;
result40 = array2table(zeros(f-s+1, iters));
for i = 1:iters
    result40(:, i) = MCMD_import(strcat(path, '\40\results\', num2str(i), '.gcmc.prp'), s, f);
end

fprintf('done!')