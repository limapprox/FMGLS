%================== Bonferroni Multivariate KPSS Tests ===================%

% <<<<< List of Inputs >>>>> %
%  T: sample size
%  n: dimensions of cross-sections
%  vu: [(nT) x 1] residuals, e.g. vuhat = vyplus-mZ*beta_fm 
%  PreMat: e.g. \widehat{\mSigma}_{\vu}^{-1}(q) 
%  'SelectStartPoint':(1) the default method construct starting points
%                          consecutively as Wagner and Hong (2016); 
%                     (2) the users can also specify it as 'COMMUTATIVE' to 
%                          follow Choi and Saikkonen (2010).

% <<<<< List of Outputs >>>>> %
%  Kstat: the Bonferroni-type KPSS test statistic, the maximum of KPSS 
%         tests constructed using subresiduals;
%  rej_rule: the p-value of Kstat times the optimal block number Mopt. Reject
%            the null rej_rule < 0.05 given the 5% sig. level; rej_rule > 1 is possible.

% <<<<< Examples of Executing this Function >>>>> %
% [Kstat,rej_rule] = multikpss_bonferroni(T,n,vuhat,invCovhat,'SelectStartPoint','CONSECUTIVE');


function [Kstat,rej_rule,bopt,Mopt] = multikpss_bonferroni(T,n,vu,PreMat,varargin)

%<<<<<<<<<<<<<<<<<<<<< Parse inputs and set defaults >>>>>>>>>>>>>>>>>>>>>%
defaultSelectStartPoint  = 'COMMUTATIVE';                 % default construction methods of starting points 
expectedSelectStartPoint = {'CONSECUTIVE','COMMUTATIVE'}; % 'CONSECUTIVE' as Wagner and Hong (2016); 'COMMUTATIVE  as Choi and Saikkonen (2010)
p = inputParser;
addParameter(p,'SelectStartPoint',defaultSelectStartPoint,@(x) any(validatestring(x,expectedSelectStartPoint))); 
parse(p,varargin{:});
SelectStartPoint = upper(p.Results.SelectStartPoint);
%=========================================================================%

%-------------------------------------------------------------------------%
%   construct Bonferroni-type KPSS Tests
bopt = MinVol(T,n,vu,PreMat,SelectStartPoint);          % select optimal block size
Mopt = floor(T/bopt);
startind = start_points(T,Mopt,bopt,SelectStartPoint);  % sequence of starting points
vk = zeros(Mopt,1);    % store values of KPSS statistics
count = 0;
for j = startind
    count     = count+1;
    vk(count) = kpss(T,n,vu,PreMat,bopt,j);
end
Kstat = max(vk);       % Bonferroni test value

%   rule of rejection
trunpot  = 100;
cdf_kpss = Fw_kpss(n,Kstat,trunpot);
if ~(cdf_kpss <= 1)
    cdf_kpss = 1;
end
rej_rule = (1-cdf_kpss)*Mopt;  % Rescale the p-values by the factor Mopt to form the 
                               % rejection rule. Reject the null if rej_rule < alpha (sig. level).

end

%-------------------------------------------------------------------------
%   select optimal block size b using minimum volatility rule, following
%   Romano and Wolf (2001), Wagner and Hong (2016).

function bopt = MinVol(T,n,vu,PreMat,SelectStartPoint)

      bmin  = floor(0.5*T^(0.5));    
      bmax  = ceil(2*T^(0.5));     
      m     = 2;                          
      Vol   = zeros(bmax-bmin+1-2*m,1);    % store the volatilities
      count = 0;
      for bi = bmin+m:bmax-m
          
          Meani  = zeros(2*m+1,1);
          Stdi   = zeros(2*m+1,1);
          counti = 0;
          for bik = bi-m:bi+m
              Mbik   = floor(T/bik);
              Sbik   = start_points(T,Mbik,bik,SelectStartPoint);  % sequence of starting points
              vkbik  = zeros(Mbik,1);
              countj = 0;
              for j = Sbik
                  countj        = countj+1;
                  vkbik(countj) = kpss(T,n,vu,PreMat,bik,j);
              end
              counti        = counti+1;
              Meani(counti) = mean(vkbik);
              Stdi(counti)  = std(vkbik);
          end
          count      = count+1;
          Vol(count) = std(Meani)+std(Stdi);          
      end
      [~,ind_min] = min(Vol);
      bopt = bmin+m+ind_min-1;
      
end


%-------------------------------------------------------------------------   
%   construct KPSS-stat
%      T: sample size
%      n: dimensions of cross-sections
%      vu: [(nT) x 1] fully modified residuals
%      PreMat: precision matrix, e.g. BIAM
%      b: block size
%      startpoint:  starting point

function kstat = kpss(T,n,vu,PreMat,b,startpoint)   
   
    j  = startpoint;
    mU = reshape(vu,[n,T]);
    mUcs = cumsum(mU(:,j:j+b-1),2);     % vector cumsum of subresiduals 
    vucs = mUcs(:);
    blk_ind = n*T-n*b+1:n*T;            % indices to choose blocks, fix the last b blocks
    kstat = b^(-2)*vucs'*PreMat(blk_ind,blk_ind)*vucs;
    
end
    
%-------------------------------------------------------------------------
%    sequence of starting points
%     T:  sample size
%     M:  number of sub-blocks
%     bT: block size
%     ConstrMethod: construction methods

function  seq_start = start_points(T,M,bT,SelectStartPoint)
   
   switch SelectStartPoint
       case 'CONSECUTIVE'
           seq_start = (1:M);
       case 'COMMUTATIVE'
           M0   = ceil(M/2);
           mInd = [2*(1:M0)'-1 2*(1:M0)'];
           seq_start = reshape([1+(mInd(:,1)-1)/2*bT T-mInd(:,2)/2*bT+1]',[1,2*M0]);
           seq_start = seq_start(1:M);      
   end  
   
end

%-------------------------------------------------------------------------
%   cdf of W_n, the limiting distribution of KPSS test for cointegration
%     n: cross-sectional dimensions
%     w: w>=0, value at which to evaluate cdf
%     trunpot: truncating point

function Fw = Fw_kpss(n,w,trunpot)
   
   vj   = 0:trunpot;                  % truncating at trunpot
   vknj = (-1).^vj.*gamma(n/2+vj)./(gamma(n/2)*gamma(vj+1));
   vlnj = n/sqrt(2)+2*sqrt(2)*vj;
   vwnj = vlnj/(2*sqrt(w));
   Fw   = 2^(n/2)*sum(vknj.*erfc(vwnj));
   
end















