% =========================================================================
% -- Optimized Coordinate Descent (ODC) Massive MU-MIMO Simulator
% -------------------------------------------------------------------------
% -- (c) 2016 Christoph Studer & Michael Wu
% -- e-mails: studer@cornell.edu and miwu@xilinx.com
% -------------------------------------------------------------------------
% -- If you use this simulator or parts of it, then you must cite our 
% -- paper: Michael Wu, Chris Dick, Joseph R. Cavallaro, Christoph
% -- Studer, "High-Throughput Data Detection for Massive MU-MIMO-OFDM 
% -- using Coordinate Descent," IEEE Transactions on Circuits and 
% -- Systems I, 2016
% =========================================================================

function OCD_sim(varargin)

% -- set up default/custom parameters

if isempty(varargin)
    
    disp('using default simulation settings and parameters...')
    
    % set default simulation parameters
    par.simName = 'ERR_128x8_16QAM'; % simulation name (used for saving results)
    par.runId = 0; % simulation ID (used to reproduce results)
    par.MR = 128; % receive antennas
    par.MT = 8; % transmit antennas (set not larger than MR!)
    par.mod = '16QAM'; % modulation type: 'BPSK','QPSK','16QAM','64QAM'
    par.trials = 1e4; % number of Monte-Carlo trials (transmissions)
    par.SNRdB_list = -10:4:10; % list of SNR [dB] values to be simulated
    par.detector = {'ZF','MMSE','SIMO','OCDBOX','OCDMMSE'}; % define detector(s) to be simulated
    par.maxiter = 2; % max number of iterations for OCD algorithms
    
else
    
    disp('use custom simulation settings and parameters...')
    par = varargin{1}; % only argument is par structure
    
end

% -- initialization

% use runId random seed (enables reproducibility)
rng(par.runId);

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ];
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i, ...
            +1-1i,+1+1i ];
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
        
end

% extract average symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.Q = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.Q,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% initialize result arrays (detector x SNR)
res.VER = zeros(length(par.detector),length(par.SNRdB_list)); % vector error rate
res.SER = zeros(length(par.detector),length(par.SNRdB_list)); % symbol error rate
res.BER = zeros(length(par.detector),length(par.SNRdB_list)); % bit error rate

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.MT,par.Q,par.trials);

% trials loop
tic
for t=1:par.trials
    
    % generate transmit symbol
    idx = bi2de(bits(:,:,t),'left-msb')+1;
    s = par.symbols(idx).';
    
    % generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1i*randn(par.MR,1));
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1i*randn(par.MR,par.MT));
    
    % transmit over noiseless channel (will be used later)
    x = H*s;
    
    % SNR loop
    for k=1:length(par.SNRdB_list)
        
        % compute noise variance (average SNR per receive antenna is: SNR=MT*Es/N0)
        N0 = par.MT*par.Es*10^(-par.SNRdB_list(k)/10);
        
        % transmit data over noisy channel
        y = x+sqrt(N0)*n;
        
        % algorithm loop
        for d=1:length(par.detector)
            
            switch (par.detector{d}) % select algorithms
                case 'ZF', % zero-forcing detection
                    [idxhat,bithat] = ZF(par,H,y);
                case 'MMSE', %  MMSE detector
                    [idxhat,bithat] = MMSE(par,H,y,N0);
                case 'SIMO', % SIMO lower bound detector
                    [idxhat,bithat] = SIMO(par,H,y,s);
                case 'OCDBOX', % Box OCD detector
                    [idxhat,bithat] = OCDBOX(par,H,y);     
                case 'OCDMMSE', % MMSE OCD detector
                    [idxhat,bithat] = OCDMMSE(par,H,y,N0);                       
                otherwise,
                    error('par.detector type not defined.')
            end
            
            % -- compute error metrics
            err = (idx~=idxhat);
            res.VER(d,k) = res.VER(d,k) + any(err);
            res.SER(d,k) = res.SER(d,k) + sum(err)/par.MT;
            res.BER(d,k) = res.BER(d,k) + sum(sum(bits(:,:,t)~=bithat))/(par.MT*par.Q);
            
        end % algorithm loop
        
    end % SNR loop
    
    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',time_elapsed*(par.trials/t-1)/60);
        tic
    end
    
end % trials loop

% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.time_elapsed = time_elapsed;

% -- save final results (par and res structure)

save([ par.simName '_' num2str(par.runId) ],'par','res');

% -- show results (generates fairly nice Matlab plot)

marker_style = {'bo-','rs--','mv-.','kp:','g*-','c>--','yx:'};
figure(1)
for d=1:length(par.detector)
    if d==1
        semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
        hold on
    else
        semilogy(par.SNRdB_list,res.BER(d,:),marker_style{d},'LineWidth',2)
    end
end
hold off
grid on
xlabel('average SNR per receive antenna [dB]','FontSize',12)
ylabel('bit error rate (BER)','FontSize',12)
axis([min(par.SNRdB_list) max(par.SNRdB_list) 1e-5 1])
legend(par.detector,'FontSize',12)
set(gca,'FontSize',12)

end

% -- set of detector functions

%% zero-forcing (ZF) detector
function [idxhat,bithat] = ZF(par,H,y)
xhat = H\y;
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% MMSE detector (MMSE)
function [idxhat,bithat] = MMSE(par,H,y,N0)
W = (H'*H+(N0/par.Es)*eye(par.MT))\(H');
xhat = W*y;
G = real(diag(W*H));
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end

%% SIMO lower bound
function [idxhat,bithat] = SIMO(par,H,y,s)
z = y-H*s;
xhat = zeros(par.MT,1);
for m=1:par.MT
    hm = H(:,m);
    yhat = z+hm*s(m,1);
    xhat(m,1) = hm'*yhat/norm(hm,2)^2;
end
[~,idxhat] = min(abs(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);
end


%% Optimized Coordinate Descent (ODC) BOX version
function [idxhat,bithat] = OCDBOX(par,H,y)

% -- initialization
[row, col] = size(H);
alpha = 0; % no regularization for BOX detector
beta = max(real(par.symbols));

% -- preprocessing
dinv = zeros(col,1);
p = zeros(col,1);
for uu=1:col
    normH2 = norm(H(:,uu),2)^2;
    dinv(uu,1) = 1/(normH2+alpha);
    p(uu,1) = dinv(uu)*normH2;
end

r = y;
zold = zeros(col,1);
znew = zeros(col,1);
deltaz = zeros(col,1);

% -- OCD loop
for iters=1:par.maxiter;
    for uu=1:col
        tmp = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
        znew(uu) = projinf(tmp,beta);
        deltaz(uu) = znew(uu)-zold(uu);
        r = r - H(:,uu)*deltaz(uu);
        zold(uu) = znew(uu);
    end
end

[~,idxhat] = min(abs(znew*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end


% project onto alpha infinity-tilde-norm ball
function sproj = projinf(s,alpha)
sr = real(s);
idxr = abs(sr)>alpha;
sr(idxr) = sign(sr(idxr))*alpha;
si = imag(s);
idxi = abs(si)>alpha;
si(idxi) = sign(si(idxi))*alpha;
sproj = sr +1i*si;
end

%% Optimized Coordinate Descent (ODC) MMSE version
function [idxhat,bithat] = OCDMMSE(par,H,y,N0)

% -- initialization
[row, col] = size(H);
alpha = 0.5*N0/par.Es; % MMSE regularization
beta = max(real(par.symbols));

% -- preprocessing
dinv = zeros(col,1);
p = zeros(col,1);
for uu=1:col
    normH2 = norm(H(:,uu),2)^2;
    dinv(uu,1) = 1/(normH2+alpha);
    p(uu,1) = dinv(uu)*normH2;
end

r = y;
zold = zeros(col,1);
znew = zeros(col,1);
deltaz = zeros(col,1);

% -- OCD loop
for iters=1:par.maxiter;
    for uu=1:col
        tmp = dinv(uu)*(H(:,uu)'*r)+p(uu)*zold(uu);
        znew(uu) = projinf(tmp,beta);
        deltaz(uu) = znew(uu)-zold(uu);
        r = r - H(:,uu)*deltaz(uu);
        zold(uu) = znew(uu);
    end
end

[~,idxhat] = min(abs(znew*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,[],2);
bithat = par.bits(idxhat,:);

end


