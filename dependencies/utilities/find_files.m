function[file_list] = find_files(pattern,root_dir,do_recursive)
%FIND_FILES returns the full path of each file matching the given pattern
%
%[file_list] = find_files(pattern,root_dir,do_recursive)
%
% Inputs
%      pattern: string that will be passed to matlab's built-in DIR function
%     root_dir: string (or cell array of strings) containing paths to
%               search [pwd]
% do_recursive: flag which determines whether to search
%               subdirectories [false]
%
% Outputs
%    file_list: cell array containing full path to every matching file
%
%Alastair Moore, March 2014


%parts of the code were inspired by thread at
% http://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab

if nargin<3 || isempty(do_recursive)
    do_recursive = 0;
end

if nargin<2 || isempty(root_dir)
    root_dir = pwd;
end

if iscell(root_dir)
    file_list = {};
    for n=1:numel(root_dir)
        %fprintf('%d - Searching root_dir: %s\n',n,root_dir{n});
        file_list = [file_list{:};find_files(pattern,root_dir{n},do_recursive)];
    end
else
    
    
    %get matching directory contents
    %fprintf('Searching root_dir: %s\n',root_dir);
    list = dir(fullfile(root_dir,pattern));
    file_list = strcat([root_dir filesep], {list.name}');
    
    if do_recursive
        %get subdirs excluding '.' and '..'
        list = dir(root_dir);
        dirnames={list([list.isdir]).name};
        dirnames=dirnames(~(strcmp('.',dirnames)|strcmp('..',dirnames)));
        
        for ii = 1:length(dirnames)
            file_list = [file_list; find_files(pattern,fullfile(root_dir,dirnames{ii}),do_recursive)];
        end
    end
end