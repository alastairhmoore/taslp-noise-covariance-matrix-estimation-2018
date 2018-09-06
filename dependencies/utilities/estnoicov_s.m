function [rvv,qq] = estnoicov_s(sig_f,rvv,qq,pp)
% ESTNOICOV_S returns the souden's method [1] to estimate the
% frequency-domain noise covariance matrix.
%
% INPUTS:
% sig_f: n_mic*n_fft frequency-domain multichannel signal vector of one frame
% rvv:	 n_mic*n_mic*n_fft noise covariance matrix
% qq:	 variables that need to be updated frame-wise
% 		 qq.ryy: signal covariance matrix
% 		 qq.spp: speech presence probability of last frame
% 
% pp:	 parameters used in [1].
%		 pp.epsi: significance level in (12)
% 		 pp.L:	 the number of periodograms for average in (13)
%		 pp.K1:	 defines the size of normalized Hann window in (19)
% 		 pp.alpha_p: smoothing parameter for multichannel speech presence probability
% 		 pp.alpha_v: smoothing parameter in (25)
%		
% OUTPUT:
% Rvv:  updated n_mic*n_mic*n_fft noise covariance matrix
% spp:	speech presence probability of one frame
%
% REFERENCES:
% [1] Souden, M., Jingdong Chen, J. Benesty, and S. Affes. an Integrated Solution for Online Multichannel Noise Tracking and Reduction.
% IEEE Transactions on Audio, Speech, and Language Processing 19, no. 7 (September 2011).
% 
% Wei Xue, Imperial College, Feb 14, 2017

%--- Check and unpack parameters ---
if nargin<4
    pp = [];
end

if isstruct(pp) || isempty(pp)
    default_fields = {'epsi','L','K1','alpha_p','alpha_v'};
    default_values = [0.01, 32, 15, 0.6, 0.92];%default values
    n_fields = length(default_fields);
	for i_field = 1:n_fields
		if ~isfield(pp,default_fields{i_field})
			eval([default_fields{i_field} '=' num2str(default_values(i_field)) ';']);
		else
			eval([default_fields{i_field} '=pp.' default_fields{i_field} ';']);
		end
	end
else
	error('please specify the parameters.');	
end

ryy = qq.ryy;
spp = qq.spp;

%--- Other default parameters ---
max_val = 1e4;
min_sap = 0.4;

%--- Allocate memory ----
[n_mic, n_fft] = size(sig_f);
phi_hat = zeros(n_fft,1);
eta_hat = zeros(n_fft,1);
beta_hat = zeros(n_fft,1);
phi_inst = zeros(n_fft,1);
sap_local = zeros(n_fft,1);

%---compute phi_0 and phi_0_h in the first frame---
if qq.phi_0 < 0
    [phi_0, phi_0_h] = solve_cum(n_mic,L,epsi);
    qq.phi_0 = phi_0;
    qq.phi_0_h = phi_0_h;
else
    phi_0 = qq.phi_0;
    phi_0_h = qq.phi_0_h;
end

w_global = hann(K1);
w_global = w_global./sum(w_global);

%--- Start Processing ----
for i_freq = 1:n_fft
	sig_f_vec = sig_f(:,i_freq);	
	ryy_i_freq = ryy(:,:,i_freq);
	rvv_i_freq = rvv(:,:,i_freq);

	rxx_i_freq = ryy_i_freq-rvv_i_freq;%2.2.a
	inv_rvv = inv(rvv_i_freq);
	phi_inst(i_freq) = real(sig_f_vec'*inv_rvv*sig_f_vec);%2.2.b
	phi_hat(i_freq) = real(trace(inv_rvv*ryy_i_freq));%2.2.c
	eta_hat(i_freq) = max(phi_hat(i_freq)-n_mic,0);% 2.2.d
	beta_hat(i_freq) = real(sig_f_vec'*inv_rvv*rxx_i_freq*inv_rvv*sig_f_vec);%2.2.e

	%---eq18---
	if phi_inst(i_freq) < phi_0 && phi_hat(i_freq) < n_mic
		sap_local(i_freq) = 1;
	elseif phi_inst(i_freq) < phi_0 && phi_hat(i_freq) >= n_mic && phi_hat(i_freq)<phi_0_h
		sap_local(i_freq) = (phi_0_h-phi_hat(i_freq))/(phi_0_h-n_mic);
	else
		sap_local(i_freq) = 0;
	end
end		

phi_global = conv(phi_inst,w_global,'same'); %eq19
phi_frame = mean(phi_inst); %eq20
sap_global = phi_global < phi_0;
sap_frame = phi_frame < phi_0;
sap = min(sap_local.*sap_global.*sap_frame,0.99);
sap = max(sap,min_sap);

%__________________________________MC-SPP and noise estimation in the first and second stage(2.4-3.2)________________________________________
spp_0 = 1./(1+(sap./(1-sap)).*(1+ eta_hat).* min(real(exp(-beta_hat./(1+eta_hat))),max_val));%using val_max to avoid inf
spp_1 = alpha_p*spp + (1-alpha_p)*spp_0;
alpha_v_hat = alpha_v + (1- alpha_v)*spp_1;


 %________________Second Stage_____________________________
for i_freq = 1:n_fft
	sig_f_vec = sig_f(:,i_freq);
	rvv_i_freq = alpha_v_hat(i_freq)*rvv(:,:,i_freq)+(1-alpha_v_hat(i_freq))*(sig_f_vec*sig_f_vec');
	ryy_i_freq = ryy(:,:,i_freq);
	
	rxx_i_freq = ryy_i_freq - rvv_i_freq; %2.2.a
	inv_rvv = inv(rvv_i_freq);
	phi_inst(i_freq) = real(sig_f_vec'*inv_rvv*sig_f_vec);%2.2.b
	phi_hat(i_freq) = real(trace(inv_rvv*ryy_i_freq));%2.2.c
	eta_hat(i_freq) = max(phi_hat(i_freq)-n_mic,0);% 2.2.d
	beta_hat(i_freq) = real(sig_f_vec'*inv_rvv*rxx_i_freq*inv_rvv*sig_f_vec);%2.2.e

	%---eq18---
	if phi_inst(i_freq) < phi_0 && phi_hat(i_freq) < n_mic
		sap_local(i_freq) = 1;
	elseif phi_inst(i_freq) < phi_0 && phi_hat(i_freq) >= n_mic && phi_hat(i_freq)<phi_0_h
		sap_local(i_freq) = (phi_0_h-phi_hat(i_freq))/(phi_0_h-n_mic);
	else
		sap_local(i_freq) = 0;
	end
end

phi_global = conv(phi_inst,w_global,'same'); %eq19
phi_frame = mean(phi_inst); %eq20
sap_global = phi_global < phi_0;
sap_frame = phi_frame < phi_0;
sap = min(sap_local.*sap_global.*sap_frame,0.99);
sap = max(sap,min_sap);

%---Final Estimation---
spp = 1./(1+(sap./(1-sap)).*(1+ eta_hat).* min(real(exp(-beta_hat./(1+eta_hat))),max_val));
alpha_v_hat = alpha_v + (1- alpha_v)*spp;

for i_freq = 1:n_fft
	rvv(:,:,i_freq) = alpha_v_hat(i_freq)*rvv(:,:,i_freq)+(1-alpha_v_hat(i_freq))*(sig_f(:,i_freq)*sig_f(:,i_freq)'); %eq24
end

%--- Re-package variables---
qq.spp = spp;

function [y0, y1] = solve_cum(n_mic,L,eps0)
F_phi = @(x) (x/L).^n_mic*L*gamma(L)/(gamma(n_mic+1)*gamma(L-n_mic+1)).*hypergeom([n_mic,L+1],n_mic+1,-x/L).*double(x>=1); %eq15
a_F = 2*n_mic*L;
B_F = (L+L-2*n_mic-1)*(L-1)/((L-2*n_mic-3)*(L-2*n_mic));
b_F = 4+(a_F+2)/(B_F-1);
% c_F = a_F*(b_F-2)/(b_F*(L-2*n_mic-1));
F_phi_h = @(x) 1.*betainc(a_F*x./(a_F*x+b_F),a_F/2,b_F/2).*double(x>=1);

y0 = solve_f(F_phi,1-eps0,1e-5);
y1 = solve_f(F_phi_h,1-eps0,1e-5);