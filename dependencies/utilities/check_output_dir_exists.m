function [] = check_output_dir_exists( dir_path )
%check_out_dir_exists creates the path if it does not exist

%check whether the path is absolute
isabsolute = false;
if isunix
    if isequal(dir_path(1),'/') || isequal(dir_path(1),'~')
        isabsolute=true;
    end
else
    %untested on windows
     if isequal(dir_path(1:2),'\\') || isequal(dir_path(2:3),':\')
        isabsolute=true;
     else
         warning('Assuming path is relative')
     end
end

if ~isabsolute
    dir_path=fullfile('.',dir_path); %prepend current directory to relative path
end

if ~exist(dir_path,'dir')
    try
    mkdir(dir_path)
    %fprintf('Created dir:\t%s\n',dir_path)
    catch err
        fprintf('Failed to create dir:\t%s\n',dir_path)
        rethrow(err)
    end
end

