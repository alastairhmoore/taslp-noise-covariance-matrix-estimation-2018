function[] = extract_mono_listening_example(filename,root_dir,outputroot,chanid)
if nargin<4||isempty(chanid)
    chanid=1;
end
[~,fn] = fileparts(filename); %remove extension
outputdir = fullfile(outputroot,fn); %append to root path
check_output_dir_exists(outputdir)
all_files = find_files(filename,root_dir,1)
for ifile = 1:length(all_files)
    [fullpath,~,~] = fileparts(all_files{ifile});
    relpath = strrep(fullpath,root_dir,''); %!! includes the leading file separator
    flatpath = strrep(relpath(2:end),filesep,'_');
    [in,fs] = audioread(all_files{ifile});
    audiowrite(fullfile(outputdir,[flatpath,'.wav']),in(:,chanid),fs)
end