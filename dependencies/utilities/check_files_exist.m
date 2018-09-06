function[bool_mat] = check_files_exist(files)
bool_mat = ones(size(files));
for ifile = 1:numel(files)
    if ~exist(files{ifile},'file')
        bool_mat(ifile) = 0;
    end
end
    