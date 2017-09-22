% function [Pano,Periodo,Amp_X,Phase_X,freq] =
function [varargout] = panorama(x,varargin)
% %
% %panorama - calculate Panorama
% %
% % Inputs:
% %    x - input signal
% %    varargin - there are 5 options 
% %        window - length of window (length of segment) - if integer,
% %                 hamming window. 
% %        noverlap - the size of overlap, must be integer between 0 and 1 
% %        nfft - Number of fft (if nfft>length(window), then zero-padding for fft) 
% %        fs - sampling freq
% %        method - 'fft': fft based panorama, 'conv': convolutional based
% panorama (ref. Autoconvolution and panorama: Augmenting second-order signal analysis - ICASSP 2014)
% %    
% %    

    %  or panorama(x,window,noverlap,nfft,method)

    
% % Outputs:
% %    Pano - averaged panorama
% %    Periodo - averaged periodogram
% %    Amp_X - summed amplitude of FFT
% %    Phase_X - phase of FFT
% %    freq - freq axis
% %
% % Other m-files required: none
% % Subfunctions: none
% % MAT-files required: none
% %
% % Written by Takashi Nakamura
% % PhD candidate in Communications and signal processing group, Imperial College London
% % takashi.nakamura14@imperial.ac.uk
% % Last revision: 10-Nov-2016
% %
% %------------- BEGIN CODE --------------

% % ================
% check for x
% % ================
% x must be column vector
if (isempty(x)||~isvector(x))
    error('panorama: arg 1 (x) must be vector.');
elseif (size(x,1)==1)
    x = x(:);
end

% number of inputs
nvarargin = length(varargin);
N = length(x); % data length

% % *******************************************
% Input variables
% % *******************************************

if nvarargin == 0
    % panorama(x)
    window = hamming(N);
    win = 'hamm';
    noverlap = 0;
    nfft = length(window);
    fs = 2;
    method = 'conv';
elseif nvarargin == 1
    % panorama(x,window) or panorama(x,method)
    
    if ischar(varargin{1})
        [window,win] = obtain_win(N);
        method = varargin{1};
    else
        [window,win] = obtain_win(varargin{1});
        method = 'conv';
    end
    
    noverlap = 0;
    nfft = length(window);
    fs = 2;
elseif nvarargin == 2
    % panorama(x,window,noverlap) or panorama(x,window,method)
    
    [window,win] = obtain_win(varargin{1});
    
    if ischar(varargin{2})
        noverlap = 0;
        method = varargin{2};
    else
        noverlap = varargin{2};
        method = 'conv';
    end
    
    nfft = length(window);
    fs = 2;
elseif nvarargin == 3
    % panorama(x,window,noverlap,nfft) or panorama(x,window,noverlap,method)
    
    [window,win] = obtain_win(varargin{1});
    noverlap = varargin{2};
    
    if ischar(varargin{3})
        nfft = length(window);
        method = varargin{3};
    else
        nfft = varargin{3};
        method = 'conv';
    end
    
    fs = 2;
    
elseif nvarargin == 4
    % panorama(x,window,noverlap,nfft,fs) or panorama(x,window,noverlap,nfft,method)
    
    [window,win] = obtain_win(varargin{1});
    noverlap = varargin{2};
    nfft = varargin{3};
    
    if ischar(varargin{4})
        fs = 2;
        method = varargin{4};
    else
        fs = varargin{4};
        method = 'conv';
    end
    
elseif nvarargin == 5
    % panorama(x,window,noverlap,nfft,fs,method)
    [window,win] = obtain_win(varargin{1});
    noverlap = varargin{2};
    nfft = varargin{3};
    fs = varargin{4};
    
    if ischar(varargin{5})
        method = varargin{5};
    else
        warning('panorama: method must be -conv-, -fft-, or -sing-');
        method = 'conv';
    end
end

% % *******************************************
% check errors
% % *******************************************

% noverlap
if length(noverlap) ~= 1 || noverlap >= 1 || noverlap < 0
    error('panorama: overlap must be integer between 0 to 1');
end

% nfft
if length(nfft) ~= 1
    warning('panorama: nfft must be integer');
    nfft = length(window);
end

% fs
if length(fs) ~= 1
    error('panorama: fs must be integer');
end

% % *******************************************
% Print parameters
% % *******************************************

% str = sprintf('N=%d,w_len=%d,win=%s,nover=%.2f,nfft=%d,fs=%d,%s',N,length(window),win,noverlap,nfft,fs,method);
% disp(str);

% % *******************************************
% Panorama
% % *******************************************

L = length(window);         % length of segment
K = round(L*noverlap);      % length of overlap
M = floor((N-L)/(L-K))+1;   % Number of segments

% need to modify freq - version for odd and even
nfft = round(nfft/2)*2;

if strcmp(method,'fft')
    % =============================
    % For DFT based panorama
    % =============================

    x_mat = zeros(L,M);
    win = window*ones(1,M);
    
    for i = 1:M
        s_index = (i-1)*(L-K) + 1;
        e_index = s_index + L - 1;
        x_mat(:,i) = x(s_index:e_index);
    end
    
    x_mat_win = x_mat.*win;
    X = fft(x_mat_win,nfft);
    X = X(1:nfft/2,:);
    
    Panorama = sum(X.^2,2)/(M*L);
    freq = 0:fs/nfft:fs/2 - fs/nfft;   
    Phase_X = atan2(imag(X), real(X));
    
elseif strcmp(method,'conv')
    % =============================
    % For convolution based panorama
    % =============================
    
    p_mat = zeros(L*2-1,M);
    
    for i = 1:M
        s_index = (i-1)*(L-K) + 1;
        e_index = s_index + L - 1;
        x_in = x(s_index:e_index).*window;
        p_mat(:,i) = conv(x_in,x_in)/L;
    end
        
    X = fft(p_mat,nfft);
    X = X(1:nfft,:);
    
    Panorama = abs(sum(X,2))/M;
    freq = 0:fs/(nfft*2):fs/2 - fs/(nfft*2);
    
elseif strcmp(method,'sing')
    % =============================
    % For Single realisation panorama
    % =============================
    
    x_mat = zeros(L,M);
    win = window*ones(1,M);
    
    for i = 1:M
        s_index = (i-1)*(L-K) + 1;
        e_index = s_index + L - 1;
        x_mat(:,i) = x(s_index:e_index);
    end
    
    x_mat_win = x_mat.*win;
    
    
    X = fft(x_mat_win,nfft);
    X = X(1:nfft/2,:);
    
    freq = 0:fs/nfft:fs/2 - fs/nfft;
    
    % obtain phase difference between previous, and following segmemt
    Phase_X = angle(X);
    diff_X_phase = [zeros(nfft/2,1),diff(Phase_X,[],2),zeros(nfft/2,1)];
    Theta = diff_X_phase(:,1:M)-diff_X_phase(:,1+1:M+1);
    Panorama = abs(sum(abs(X).^2 .*cos(Theta),2))/(L*M);
    
else
    error('panorama: method must be -conv-, -fft-, or -sing-');
end

varargout{1} = Panorama;

if nargout == 2
    varargout{2} = freq;
elseif nargout == 3
    varargout{2} = freq;
    if strcmp(method,'conv')
        warning('panorama: conv method only have output 1 and 2');
        varargout{3} = [];
    else
        varargout{3} = Phase_X;
    end
elseif nargout == 4
    varargout{2} = freq;
    if strcmp(method,'conv')
        warning('panorama: conv method only have output 1 and 2');
        varargout{3} = [];
        varargout{4} = [];
    elseif strcmp(method,'fft')
        warning('panorama: fft method only have output 1-3');
        varargout{3} = Phase_X;
        varargout{4} = [];
    else
        varargout{4} = Theta;
    end
end

end

% obtain window
function [window,win] = obtain_win(input)
if length(input) == 1
    window = hamming(input);
    win = 'hamm';
else
    if (size(input,1)==1)
        window = input(:); % window must be column vector
    end
    win = 'option';
end
end

% %------------- END OF CODE --------------
