function[out_path] = writewav_with_gain_factors(z,fs,output_dir,tags)

% forces writewav to scale the output to be in range +/-1 and save as 32
% bit floating point values. Use readwave with 'g' option to recover the absolute values
mode_str = 'gesvL'; 

nSigs = length(z);
out_path = cell(nSigs,1);
if size(z)~=size(tags)
    error('There must be the same number of tags as there are signals')
end

check_output_dir_exists(output_dir);

try
    for n=1:length(z)
        out_path{n} = fullfile(output_dir,sprintf('%d_%s.wav',n,tags{n}));
        writewav(z{n},fs,out_path{n},mode_str);
    end
    fprintf('\n\nProcessed audio saved to %s\n',output_dir);
catch
    error('Error saving files to %s',output_dir);
end

