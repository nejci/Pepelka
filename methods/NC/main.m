function ADDED_PATH=main()
% main;
% 
% Initializes matlab paths to subfolders
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.

files = dir(cd);

ADDED_PATH=[];

for i=1:length(files)
    if files(i).isdir && strcmp(files(i).name,'.') == 0  && strcmp(files(i).name,'..') == 0
        new_path=[cd '/' files(i).name];
		addpath(new_path);
		ADDED_PATH=[ADDED_PATH, new_path,';'];
    end
end
