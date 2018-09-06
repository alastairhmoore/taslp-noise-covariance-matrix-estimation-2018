function[w,issingular] = fcn_20170222_01_mvdr_weights(R,d,debug_level)
%R is spatial covariance of stft of received signal [nChans,nChans,nFreq,nFrames]
%d is steering vector for each TF bin [nChans,nFreq,nFrames]
%debug_level controls how to handle singular or ill-conditioned matrices
% 0: ignore
% 1: print warning to console (as standard)
% 2: catch and allow debugging
% 3: throw error

% define the msg id for singular matrix inversion
MSGID_NEAR_SINGULAR = 'MATLAB:nearlySingularMatrix';
MSGID_SINGULAR = 'MATLAB:singularMatrix';
if 0
  % handy to have these lines for debugging
  warning('on', MSGID_NEAR_SINGULAR);
  warning('on', MSGID_SINGULAR);
end


if nargin<3 || isempty(debug_level)
    debug_level = 1;
end

if nargout==1 && debug_level==0
    % don't care about singular values - suppress warnings
    prev_warn_state(1) = warning('off', MSGID_NEAR_SINGULAR);
    prev_warn_state(2) = warning('off', MSGID_SINGULAR);
end
if nargout == 1 && debug_level==1
    % everything can stay as it is, but define the state to be sure
    prev_warn_state(1) = warning('on', MSGID_NEAR_SINGULAR);
    prev_warn_state(2) = warning('on', MSGID_SINGULAR);
end
% any other states and we turn warnings into errors so they can be caught and
% processed accordingly
% N.B. based on neat trick documented at http://undocumentedmatlab.com/blog/trapping-warnings-efficiently
if ~exist('prev_warn_state','var')
    prev_warn_state(1) = warning('error', MSGID_NEAR_SINGULAR);
    prev_warn_state(2) = warning('error', MSGID_SINGULAR);
end
% ensure warnings are returned to their previous states on exiting the
% function (whether by error or completion)
%c = onCleanup(@()eval('warning(prev_warn_state)'));
c = onCleanup(@()cleanup(prev_warn_state));

% validate sizes
[nChansR1,nChansR2,nFreqR,nFramesR] = size(R);
[nChansd,nFreqd,nFramesd] = size(d);

nChans = unique([nChansR1,nChansR2,nChansd]);
nFreq = unique([nFreqR,nFreqd]);
nFrames = unique([nFramesR,nFramesd]);
if numel(nChans)~=1 || numel(nFreq)~=1 || numel(nFrames)~=1
    error('Dimensions of R and d are not consistent')
end

%preallocate
w = zeros(nChans,nFreq,nFrames);
issingular = zeros(nFreq,nFrames);








%for clarity and simplicity, loop over time and frequency indices
for ifreq = 1:nFreq
    for iframe = 1:nFrames
        try
            % ---------
            invRd = R(:,:,ifreq,iframe) \ d(:,ifreq,iframe);
            % ---------
        catch ME
            % handle inversion of singular matrix
            if strcmp(ME.identifier, MSGID_NEAR_SINGULAR) || ...
                    strcmp(ME.identifier, MSGID_SINGULAR)
                
                issingular(ifreq,iframe) = 1;
                
                switch debug_level
                    case 0
                        % do nothing
                    case 1
                        % print the warning
                        warning(ME.message)
                    case 2
                        % print the warning and stop execution
                        warning(ME.message)
                        keyboard;
                    case 3
                        rethrow(ME);
                end
                
                % if we get to here then we want to continue, but the matrix
                % inverse errored out so need to actually compute it
                warning('off', MSGID_NEAR_SINGULAR);
                warning('off', MSGID_SINGULAR);
                
                invRd = R(:,:,ifreq,iframe) \ d(:,ifreq,iframe);
                
                warning('error', MSGID_NEAR_SINGULAR);
                warning('error', MSGID_SINGULAR);
                
            else
                rethrow(ME);
            end
        end
        
        % ---------
        w(:,ifreq,iframe) = invRd / (d(:,ifreq,iframe)' * invRd);
        % ---------
    end
end

function[] = cleanup(prev_warn_state)
% the warnings
warning(prev_warn_state)
