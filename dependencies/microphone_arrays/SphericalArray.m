classdef (HandleCompatible) SphericalArray
    
    
    
end


 function[Hnm,f] = getFreqDomainArrayManifoldAsComplexSH(obj,Nharm,n_fft)
                    switch obj.nativeFormat
                        case 'analytical'
                            % for analytical array we can sample the manifold at arbitrary resolution
                            % although there are fast methods for computing the SH trasnform,
                            % here we do a brute force approach for each
                            
                            nSH = (Nharm+1)^2;
                            nDirReq = nSH; %minimum!
                            
                            if nDirReq>794
                                error('Revisit this')
                            end
                            load sph_dist_794.mat
                            src_az = az_794;
                            src_inc = inc_794;
                            
                            Y = sphBasis(src_az,src_inc,Nharm);
                            Yinv = pinv(Y);
                            
                            [ir,t0] = obj.getIRForSrcDoa(src_az,src_inc);
                            
                            if nargin<3 || isempty(n_fft)
                                n_fft = size(ir,1);
                            end
                            
                            if size(ir,1)>n_fft
                                error('Need to handle this')
                            end
                            H = rfft(ir,n_fft,1);
                            [nFreq,nMic,nSrc] = size(H);
                            %need 1st dimension to be wrt src_az
                            H = permute(H,[3 1 2]); %[nSrc,nFreq,nMic]
                            Hnm = zeros(nFreq,nSH,nMic);
                            for imic = 1:nMic
                                for ifreq = 1:nFreq
                                    Hnm(ifreq,:,imic) = Yinv * H(:,ifreq,imic);
                                end
                            end
                            f = (0:fix(n_fft/2)).' * (obj.fs/n_fft);
                            
                            %% debug
                            if 1
                                % use SH coefficients to reconstrunct directional
                                % response
                                % Y is [nSrc,nSH]
                                % Hnm is [nFreq,nSH,nMic]
                                % want H_hat as [nSrc,nFreq,nMic] for compatibility
                                % with the above
                                % move nSH into 4th dimension and sum across it
                                H_hat = sum(bsxfun(@times,...
                                    permute(Y,[1 3 4 2]),...            % [nSrc 1     1    nSH]
                                    permute(Hnm,[4 1 3 2])),...         % [1    nFreq nMic nSH]
                                    4);
                                comp = {'mag','phase'}
                                quant = {'orig','recon'}
                                fsel = [100 1000 7500]
                                lfs = length(fsel);
                                
                                for ic = 1:length(comp)
                                    for iq = 1:length(quant)
                                        figure
                                        setFigureSize(gcf,[30 20],'centimeters');
                                        for ii = 1:lfs
                                            [thisf,iif] = find_nearest(fsel(ii),f);
                                            for imic = 1:nMic
                                                ph(ii,imic) = subplot(lfs,nMic,(ii-1)*nMic + imic);
                                                switch quant{iq}
                                                    case 'orig'
                                                        D = H(:,iif,imic);
                                                    case 'recon'
                                                        D = H_hat(:,iif,imic);
                                                end
                                                switch comp{ic}
                                                    case 'mag'
                                                        D = 20*log10(abs(D));
                                                    case 'phase'
                                                        D = angle(D);
                                                end
                                                scatter(radtodeg(src_az),radtodeg(src_inc),20,D,'filled');
                                                title(sprintf('%2.2f Hz',thisf));
                                                xlabel('Azimuth [deg]')
                                                ylabel('Inclination [deg]')
                                            end
                                        end
                                        setappdata(gcf,'lp',linkprop(ph,{'xlim','ylim','clim'}));
                                        switch comp{ic}
                                            case 'mag'
                                                colormap(v_colormap('v_thermliny'))
                                            case 'phase'
                                                phasemap
                                        end
                                        ahm_print_to_pdf(gcf,sprintf('Array manifold %s %s',quant{iq},comp{ic}));
                                    end
                                end
                                
                                %keyboard
                            end
                        otherwise
                            error('Function not yet implemented for %s arrays',obj.nativeFormat);
                    end