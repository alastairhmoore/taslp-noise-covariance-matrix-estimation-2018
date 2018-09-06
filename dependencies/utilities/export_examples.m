function[] = export_examples(in_dir,out_dir)

in_files = find_files('*.wav',in_dir);
if isempty(in_files)
    error('No files found')
end
check_output_dir_exists(out_dir);

for ifile=1:length(in_files)
    [~,fn,ext] = fileparts(in_files{ifile});
    [x,fs] = audioread(in_files{ifile});
    audiowrite(fullfile(out_dir,[fn ext]),...
        x(:,1:min(size(x,2),2)),fs);
end