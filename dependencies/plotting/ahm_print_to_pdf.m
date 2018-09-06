function[] = ahm_print_to_pdf(fh,filename)

[fpath,fname,ext] = fileparts(filename);
if ~isempty(fpath)
    check_output_dir_exists(fpath);
end
if isempty(ext)
    ext = '.pdf';
end

%% size manipulations
%fillFigure;
px = 0.01; %fudge factors to avoid clipping, set to zero to restore old behaviour
py = 0.01;
set(fh,'Units','centimeters')
set(fh,'PaperUnits','centimeters');
pos = get(fh,'position');
x=pos(3);
y=pos(4);
set(fh,'PaperPosition',[px*x, py*y, x, y]);
set(fh,'PaperSize',[(1+2*px) * x, (1+2*px) * y]);

%% output
filepath = fullfile(fpath,[fname ext]);
switch ext
    case '.pdf'
        % extra bit to embed the fonts
        switch computer
            case 'MACI64'
                if 0
                    tmpfile = [tempname ext];
                    print(fh,tmpfile,'-dpdf');
                    exe_str = '/opt/local/bin/gs';
                    %exe_str = '/usr/local/bin/gs';
                    template_str = ['%s -dSAFER -dNOPLATFONTS -dNOPAUSE -dBATCH -sDEVICE=pdfwrite ', ...
                        ' -dCompatibilityLevel=1.5 -dPDFSETTINGS=/printer ', ...
                        '-dMaxSubsetPct=100 -dSubsetFonts=true -dEmbedAllFonts=true -sOutputFile=%s -f %s'];
                    
                    [st,res] = system(sprintf(template_str,exe_str,filepath,tmpfile));
                    if st
                        disp(res)
                    end
                else
                    print(fh,filepath,'-dpdf');
                end
            otherwise
                print(fh,filepath,'-dpdf');
                
        end
    case '.png'
        print(fh,filepath,'-dpng','-r300');
end