%%
clear; clc; close all
opts = delimitedTextImportOptions("NumVariables", 15);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = [" ", ":"];

% Specify column names and types
opts.VariableNames = ["ITEM", "VarName2", "TIMESTEP", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "Var13", "Var14", "Var15"];
opts.SelectedVariableNames = ["ITEM", "VarName2", "TIMESTEP", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12"];
opts.VariableTypes = ["double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string"];

% Specify file level properties
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var13", "Var14", "Var15"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["TIMESTEP", "VarName4", "VarName5", "VarName6", "VarName7", "Var13", "Var14", "Var15"], "EmptyFieldRule", "auto");

% Import the data
ch41 = readtable("C:\Users\btapi\Downloads\ch4_1.txt", opts);
ch42 = readtable("C:\Users\btapi\Downloads\ch4_2.txt", opts);
ch43 = readtable("C:\Users\btapi\Downloads\ch4_3.txt", opts);
ch44 = readtable("C:\Users\btapi\Downloads\ch4_4.txt", opts);
ch4 = readtable("C:\Users\btapi\Downloads\ch4.txt", opts);


clear opts
%%
x1 = table2array(ch41(:,4));
y1 = table2array(ch41(:,5));
z1 = table2array(ch41(:,6));
x2 = table2array(ch42(:,4));
y2 = table2array(ch42(:,5));
z2 = table2array(ch42(:,6));
x3 = table2array(ch43(:,4));
y3 = table2array(ch43(:,5));
z3 = table2array(ch43(:,6));
x4 = table2array(ch44(:,4));
y4 = table2array(ch44(:,5));
z4 = table2array(ch44(:,6));
init1 = [x1(1),y1(1),z1(1)];
init2 = [x2(1),y2(1),z2(1)];
init3 = [x3(1),y3(1),z3(1)];
init4 = [x4(1),y4(1),z4(1)];
disp1 = ((x1-init1(1)).^2+(y1-init1(2)).^2+(z1-init1(3)).^2).^(1/2);
disp2 = ((x2-init2(1)).^2+(y2-init2(2)).^2+(z2-init2(3)).^2).^(1/2);
disp3 = ((x3-init3(1)).^2+(y3-init3(2)).^2+(z3-init3(3)).^2).^(1/2);
disp4 = ((x4-init4(1)).^2+(y4-init4(2)).^2+(z4-init4(3)).^2).^(1/2);
msd1 = disp1.^2; 
msd2 = disp2.^2;
msd3 = disp3.^2;
msd4 = disp4.^2;
msd = (msd1+msd2+msd3+msd4)/4;

spacing = 100;
len = length(msd)*spacing;

axis = 100:spacing:len;
%%
figure
plot3(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
figure
scatter(axis,[msd1,msd2,msd3,msd4,msd])
figure
loglog(axis,[msd1,msd2,msd3,msd4,msd])
delta = 1;
tau = 1;
for c = 1:1
    
end