function c = run_asymp_pdc(X, varargin)

defaults = struct();
defaults.alg = 1;
defaults.criterion = 1;
defaults.maxIP = 30;
defaults.VARadequacy_signif = 0.05;
defaults.metric = 'info';
defaults.alpha = 0;
defaults.fs = 1e3;
defaults.nFreqs = 128;
defaults.verbose = 0;
defaults.dir_method = 'pdc';

params = shared_utils.general.parsestruct( defaults, varargin );

dir_method = validatestring( params.dir_method, {'pdc', 'dtf'}, mfilename, 'dir_method' );

%        metric   'euc'  - Euclidean     -> original PDC or DTF;
%                 'diag' - diagonal      -> gPDC or DC;
%                 'info' - information   -> iPDC or DTF;

alg = params.alg; %<***> Nuttall-Strand (alg=1) algorithm, our experience indicates that
                  %      N-S is a good and robust MAR model estimation algorithm.

criterion = params.criterion; %<***> AIC, Akaike information criterion (Our preferred one 
                              %      along with BIC)
  
maxIP = params.maxIP; % maxIP - externally defined maximum IP %<***>

% Adequacy of MAR model estimate
VARadequacy_signif = params.VARadequacy_signif; % VAR model adequacy significance level
% Selection of metric used in PDC or DTF estimate

metric = params.metric;

% Significance level for PDC2 null hypothesis testing, GCT and iGCT 
alpha = params.alpha;        %<***> Significance level for PDC2 null hypothesis testing,
                     %
                     % Note: if alpha == 0, no asymptotic statistics 
                     % calculation is performed and ASYMP_PDC (see bellow) 
                     % will only returns PDC. This option is interesting 
                     % if you want faster PDC calculation.

fs = params.fs;
nFreqs = params.nFreqs;

w = fs*(0:(nFreqs-1))/2/nFreqs; % frequency scale

verbose = params.verbose;

%%

[IP,pf,A,pb,B,ef,eb,vaic,Vaicv] = mvar(X,maxIP,alg,criterion, verbose);

if ( strcmp(dir_method, 'pdc') )
  c = asymp_pdc(X,A,pf,nFreqs,metric,alpha);
  
elseif ( strcmp(dir_method, 'dtf') )
  c = asymp_dtf(X,A,pf,nFreqs,metric,alpha); % for DTF analysis
  c.pdc = c.dtf;
  c.pdc2 = c.dtf2;
  
else
  error( 'Unhandled dir_method "%s".', dir_method );
end

c.f = w;

end