function [PDC_Out, PDC_Avg, f] = partialDirectedCOH(LFPData)
%PARTIALDIRECTEDCOHFUNC Summary of this function goes here
%   Input is matrix LFPData with the size of nTrials x nSamples x nChans
%   Where nTrials is the number of trials, nSamples is the length of each
%   trial in samples and nChans is the number of channels
% maybe change to parfor?
%% LFP data Settings
[nTrials, nSamples, nChans] = size(LFPData);
Fs = 1000; % Sampling frequency
Fmax = 90;        % Maximum frequency band in the T-F GOPDC measure
nFreqBins = 2*Fmax;      % Number of frequency bins from 0 to Fmax
f = linspace( 0, Fmax, nFreqBins );

PDC_Out = zeros(nChans,nChans,nFreqBins,nSamples,nTrials);

%% Loop through trials
for t = 1 : nTrials
    
    y = squeeze(LFPData(t,:,:)); % each epoch
    
    %% DEKF for time-varying MVAR model estimation
    inp_model.data = y;
    inp_model.order = 5;             % Model order
    [A,C] = DEKF3(inp_model);        % Estimated time-varying parameters, A = [A1 A2 ... Ar] --> A: CH x CH*p x L
    
    %% Compute connectivity measures including GOPDC
    [GPDC,OPDC,PDC,GOPDC,S] = PDCFunc(A,C,inp_model.order,Fs,Fmax,nFreqBins);
    
    %% Save the PDC into matrix
    PDC_Out(:,:,:,:,t) = GOPDC(:,:,:,:);
    
    
end

PDC_Avg = mean(PDC_Out,5);

end


function [GPDC,OPDC,PDC,GOPDC,S] = PDC(A,C,p_opt,Fs,Fmax,Nf)
% A = [A1 A2 ... Ap]
% Nf: Number of frequency points
% Fmax: maximum frequency limit (should be less than Fs/2)

[CH,L,T] = size(A);

PDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence
OPDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence
GPDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence
GOPDC   = zeros(CH,CH,Nf,T); % Orthogonalized Partial Directed Coherence

S   = zeros(CH,CH,Nf,T); % Partial Directed Coherence
% f = (0:Nf-1)*((Fs/2)/Nf); % Frequency span ---> original code
f = (0:Nf-1)*(Fmax/Nf); % Frequency span
z = 1i*2*pi/Fs;

for t = 1 : T % Number of time points
    clear A2
    A2 = [eye(CH) -squeeze(A(:,:,t))];
    
    C_tmp = squeeze(C(:,:,t));
    Cd = diag(diag(C_tmp)); % Cd is useful for calculation of DC
    invCd = abs(inv(Cd));% invCd is useful for calculation of GPDC
    if(sum(sum(isinf(invCd)))>0)
        invCd = zeros(CH);
    end
    
    for n = 1 : Nf
        A_f = zeros(CH);
        for k = 1 : p_opt+1
            A_f = A_f + (A2(:,k*CH+(1-CH:0))*exp(z*(k-1)*f(n))); % sum(Ar(f)), r=0:p --> A(f) (if B==I)
        end
        
        H_f = inv(A_f);
        
        %% Power Spectral Density matrix
        S(:,:,n,t) = (H_f * C_tmp * H_f');
        
        for ii = 1:CH
            a_ii = squeeze(A_f(:,ii)); % ii'th column of A(f) --> ak
            a_denum(ii) = sqrt(a_ii'*a_ii); % hermitian(ak) * ak (k=1:CH)
            a2_denum(ii) = sqrt(a_ii'*invCd*a_ii); % for the GPDC - uses diagonal covariance information
        end
        
        PDC(:,:,n,t)  = A_f./a_denum(ones(1,CH),:); % Aij/(sqrt(aj'*aj))
        OPDC(:,:,n,t)  = (real(A_f).*(imag(A_f))./(a_denum(ones(1,CH),:)).^2); % Aij/(sqrt(aj'*aj))
        
        GPDC(:,:,n,t)  = sqrt(invCd)*A_f./a2_denum(ones(1,CH),:); % Aij/(sqrt(aj'*aj))
        GOPDC(:,:,n,t)  = (invCd)*(real(A_f).*(imag(A_f))./(a2_denum(ones(1,CH),:)).^2); % Aij/(sqrt(aj'*aj))
        
    end
    % t
end


end

function [A,Sigma_e] = DEKF3(inp_model)
% A: Estimated time-varying parameters, A = [A1 A2 ... Ar]
% inp_model.data: (L x CH) data matrix
% inp_model.order: Model order
% inp_model.param{i}.vector: Data vector of the i'th time-varying parameter
% inp_model.param{i}.name: Name of the i'th time-varying parameter
% UC: Update Coefficient or Forgetting Factor Coefficient

%% Written by: Amir Omidvarnia
y = inp_model.data;
p = inp_model.order;
L = size(y,1); % Number of samples

%%
y = y';
M = size(y,1);                    % Number of states (here, M = N)
LEN = size(y,2);                  % Number of the multivariate observations

%% Initial parameters for Dual Extended Kalman Filter
%%%%% (EKF 1)
xh = zeros(M*p,LEN);            % (EKF 1) Initial a-posteriori states (Mp x 1)
Px = .1*eye(M*p);               % (EKF 1) Initial a-posteriori state covariance matrix
% R = zeros(M,M,LEN);             % (EKF 1,2) Measurement error covariance matrix
% R(:,:,p+1) = 100*eye(M);          % (EKF 1,2) Initial observation noise covariance matrix
R = eye(M);
Q = 10*eye(M);                  % (EKF 1,2) Initial process noise covariance matrix
B = zeros(M*p,M);               % (EKF 1) Relationship between states and the process noise ( x(k) = F[x(k-1)] + B*v(k) )
B(1:M,:) = eye(M);              % (EKF 1) B = [I 0 ... 0]'
C = B';                         % (EKF 1,2) Measurement matrix (identity matrix, C = B')
% Kx = zeros(M*p,M,LEN);

% ah  = zeros(M*M*p,LEN);         % (EKF 2) Initial a-posteriori parameters estimaes
Ah = zeros(M*p,M*p,LEN);        % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
Ah(1:M,1:M*p,p) = .5*randn(M,M*p);     % (EKF 2) Initial a-posteriori parameters estimates (matrix form of 'ah' plus identity matrices)
for r = 2 : p
    Ah((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M, p) = eye(M);
end

%%%%% EKF 2
Pa = eye(M*M*p);                % (EKF 2) Initial a-posteriori parameters covariance matrix
% Ka = zeros(M*M*p,M,LEN);

for r = 1 : p
    xh((r-1)*M+1:r*M,p+1) = y(:,p-r+1);
end

Sigma_e = zeros(M,M,LEN);
%% DEKF starts ....
for i = p+1 : LEN
    
    [J_x J_A] = MVAR_JacCSD(Ah(:,:,i-1),xh(:,i-1),p); % xh(k) = F(A(k-1) * xh(k-1)) = Ah(k-1) * xh(k-1)
    Ah_ = Ah(:,:,i-1);                                 % Ah_(k) = Ah(k-1)
    
    %% EKF 1 ---> States estimation
    %---------- Time Update (EKF1) ----------
    Rv = B * Q * B';                          % According to Haykin's book
    xh_ = Ah_ * xh(:,i-1);                    % xh_(k) = A_h(k-1) * xh(k-1)
    Px_ = J_x * Px * J_x' + Rv;               % Px_(k) = A_h(k-1) * Px(k-1) * A_h(k-1)' + B * Q * B'
    
    %---------- Measurement Update (EKF1) ----------
    Rn = R; %R(:,:,i-1);                                            % According to Haykin's book
    Kx = Px_ * C' * inv(C * Px_ * C' + Rn);            % Kx(k)  = Px_(k) * C' * inv(C * Px_(k) * C' + R)
    Px = (eye(M*p) - Kx * C) * Px_;                    % Px(k)  = (I - Kx(k) * C) * Px_(k)
    e = y(:,i) - C * xh_;          % inov(k) = y(k) - C * Ah_(k) * xh(k-1)
    xh(:,i) = xh_ + Kx * e;           % xh(k)  = xh_(k) + Kx(k) * (y(k) - C * xh_(k))
    
    
    %% EKF 2 ---> Parameters estimation
    %---------- Time Update (EKF2) ----------
    ah_ = reshape(Ah_(1:M,:)',M*M*p,1);                % ah_ = vec(Ah(k-1))
    Rr = .02*Pa;                                      % Rr = lambda * Pa(k-1)
    Pa_ = Pa + Rr;                                     % Pa_(k) = Pa(k-1) + Rr
    %---------- Measurement Update (EKF2) ----------
    %%%%%% Compute DfDa and H (H(k) = C*DfDa(k-1))
    H = C * (-J_A);                                    % J_A = -DfDa; ---> Why?
    %%%%%%%%%%%%%%%
    Re = (Rn + Q);                                        % According to Haykin's book
    Ka = Pa_ * H' * inv(H * Pa_ * H' + Re);            % Ka(k) = Pa_(k) * H(k) * inv(H(k) * Pa_(k) * H(k)' +R + Q)
    Pa = (eye(M*M*p) - Ka * H) * Pa_;                  % Pa(k) = (I - Ka(k) * H(k)) * Pa_(k)
    ah = ah_ + Ka * (y(:,i)- C * Ah_ * xh(:,i-1));                    % ah(k) = ah_(k) + Ka(k) * (y(k) - yh_(k))
    
    % Re-arrange vector ah(k) into the matrix Ah(k)
    Ah(1:M,1:M*p,i) = reshape(ah,M*p,M)';
    for r = 2 : p
        Ah((r-1)*M+1:r*M,(r-2)*M+1:(r-1)*M, i) = eye(M);
    end
    
    %     Update the estimation error covariance matrix for dDTF, Partial
    %     Coherence and adaptive AIC
    Sigma_e(:,:,i) = C * Px_ * C' + R; % Ref: Advanced digital Signal Processing and noise redution (4th ed.), Saeed Vaseghi, Eq. 7-18
    %       i
end

A = Ah(1:M,:,:); % Estimated time-varying parameters, A = [A1 A2 ... Ar]
end

%% How to plot the results

% %%%%%
% figure, % ---> GOPDC
% s1 = 0;
% clear mask
% for i = 1 : nChans
%     for j = 1 : nChans
%         s1 = s1 + 1;
%         h = subplot(nChans,nChans,s1);
%         set(h,'FontSize',15,'FontWeight','bold');
%
%         PDC_tmp = abs(PDC_Avg(i,j,:,:));
%
%         img = squeeze(PDC_tmp);
%         if(i==j)
%             img = zeros(size(PDC_Avg,3),size(PDC_Avg,4));
%         end
%         imagesc([0 400],[2 Fmax],img)
%         caxis([0 .05])
%         set(h,'YDir','normal')
%         grid on
%
%         if(i==nChans && j==ceil(nChans/2))
%             xlabel('Time (msec)','Fontsize',20)
%         end
%         if(i==ceil(nChans/2) && j==1)
%             ylabel('Frequency (Hz)','Fontsize',20)
%         end
%
%     end
% end
% h2 = colorbar;
% set(h2, 'Position', [.92 .11 .03 .8150],'FontSize',20,'FontWeight','bold')

