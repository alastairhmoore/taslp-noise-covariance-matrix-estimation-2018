function[scale_factor,out_path] = audiowrite_in_subfolders_with_joint_normalisation(z,fs,abs_max_val,output_dir,basename,tags)

nSigs = length(z);
out_path = cell(nSigs,1);
if size(z)~=size(tags)
    error('There must be the same number of tags as there are signals')
end
for n=1:length(z)
    check_output_dir_exists(fullfile(output_dir,sprintf('%d_%s',n,tags{n})));
end

% find scale factor to do joint normalisation
extreme_vals = zeros(numel(z),2);
for n=1:numel(z)
    extreme_vals(n,:) = prctile(z{n}(:),[0 100]);
end
[~,scale_factor] = normalise(extreme_vals(:)./abs_max_val,16);

try
    for n=1:length(z)
        %wavwrite(0.8*normalise(z{n},16),fs,16,fullfile(output_dir,sprintf('%s_%d_%s.wav',basename,n,tags{n})))
        out_path{n} = fullfile(output_dir,sprintf('%d_%s',n,tags{n}),sprintf('%s.wav',basename));
        audiowrite(out_path{n},...
            scale_factor*z{n},fs,...
            'BitsPerSample',16);
    end
    fprintf('\n\nProcessed audio saved to %s\n',output_dir);
    % write a text file to keep track of the scale factor
    sfdir = fullfile(output_dir,'scale_factors');
    check_output_dir_exists(sfdir)
    fid = fopen(fullfile(sfdir,sprintf('%s.txt',basename)),'w');
    fprintf(fid,'%1.10f',scale_factor);
    fclose(fid);
catch
    error('Error saving files to %s',output_dir);
end

