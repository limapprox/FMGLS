%############ HAC: Andrews (1991)  ############%

% This function gives a long-run covariance estimate for applying FMOLS.

% <<<<< List of Inputs >>>>> %
% *1.  Xi = [vxi_1'
%            vxi_2'
%             ...
%            vxi_T'] 
%      is a [T x (m+n)] matrix with vxi_t = [vv_t';vu_t']', where vv_t is 
%      a [m x 1] vector, and vu_t is a [n x 1] vector;
% *2. 'Kernel': 'TR','BT','PZ','TH','QS';
% *3. 'Bandwidth': (1) the default bandwidth selection employs data-dependent
%                 rule of Andrews (1991) using OLS regression for AR(1)
%                 model ('AR1OLS'). The user is free to choose the other 
%                 methods as well including 'AR1','AR1MLE','ARMA11';
%                 (2) the user is also free to set the bandwidth by using:
%                 ['Bandwidth', 0.5], say;
% *4. 'Prewhiten': (1) the default programme does not apply prewhitening
%                 procedure which is developed by Andrews and Monahan (1992);
%                  (2) if the user would like to apply this procedure, 
%                 specify ['Prewhiten',p] for some p>0. Then it can
%                 prewhiten the process by a VAR(p) process.

% <<<<< List of Outputs >>>>> %
% 
% LongRunCov : \widehat{\mOmega}_{xi}, the long-run vairance of [\vu_t,\vvt]' 
% OneSidedLRV: \widehat{\mDelta}_{xi}, the one-sided long-run vairance of [\vu_t,\vvt]'
% BiasCorr_SimulatedEst: {\widehat{\mSigma}-\widehat{\mDelta}'}_{vu}, for bias correction in NLS

% <<<<< Examples of Executing this Function >>>>> %
% The following examples are based on Wagner & Hong (2016, page 1303). It
% applies Bartlett ('BT') or Quadradic Spectral ('QS') kernel. The bands
% are selected by the data-dependent procedure of Andrews (1991). There is
% no prewhitening procedure of Andrews and Monahan (1992) implemented.
%    Ex1: LongRunCov = AndHAC(Xi,'Kernel','BT');
%    Ex2: LongRunCov = AndHAC(Xi,'Kernel','QS');
% To obtain the one-sided long-run covariance, try
%    Ex3: [LongRunCov,OneSidedLRV] = AndHAC(Xi,'Kernel','BT');

% If the user would like to determine the bandwidth say 0.5, try
%    Ex4: LongRunCov = AndHAC(Xi,'Kernel','BT','Bandwidth',0.5);
% or using the data-dependent rule but using other methods say MLE, try
%    Ex5: LongRunCov = AndHAC(Xi,'Kernel','BT','Bandwidth','AR1MLE');

% If the user would like to prewhiten the process by a VAR(1) process, try
%    Ex6: LongRunCov = AndHAC(Xi,'Kernel','BT','Prewhiten',1);



function [LongRunCov,OneSidedLRV,BiasCorr_SimulatedEst] = AndHAC(Xi,varargin)

   [T,N]=size(Xi);      % T is sample size, N:= m+n


%<<<<< Parse inputs and set defaults >>>>>%

   defaultKernel     = 'BT';      % default Kernel: Bartlett kernel
   expectedKernel    = {'TR','BT','PZ','TH','QS'}; 
   defaultPrewhiten  = 0;         % 0: non-prewhitening; p>0: prewhitened by VAR(p)
   
   
   p = inputParser;
   addParameter(p,'Kernel',defaultKernel,@(x) any(validatestring(x,expectedKernel))); 
   addParameter(p,'Bandwidth',[],@CheckBandwidth);
   addParameter(p,'Prewhiten',defaultPrewhiten,@(x) isnumeric(x) && x>=0 && mod(x,1) == 0); 
   parse(p,varargin{:});
   
   Kernel = upper(p.Results.Kernel);
   b = p.Results.Bandwidth;
   p = p.Results.Prewhiten;   
   
   
%-------------------------------------------------------------------------
%  Prewhitening procedure based on Andrews and Monahan (1992)

   if p > 0            % prewhitened by VAR(p) process
       if p > T-1
           error(message('the lag order in prewhitening procedure cannot be larger than T-1'))
       end
       
       VARp = varm('AR',nancell(N,p),'Constant',zeros(N,1));
       Vx   = Xi(1:p,:);
       Vy   = Xi(p+1:end,:);
       [VARfit,~,~,V] = estimate(VARp,Vy,'Y0',Vx);
       VARcoeffs      = VARfit.AR;
   elseif p == 0       % no prewhitening
       V = Xi;
   end
  
%-------------------------------------------------------------------------
% Get an optimal bandwidth if it is unspecified by the user
%  Selection method includes: 'AR1','AR1OLS','AR1MLE','ARMA11'

   if isempty(b) || ischar(b)
       if isempty(b)
           model = 'AR1OLS';            % default method
       elseif ischar(b)
           model = b;  
       end
       try
           b = getBandwidth(V,Kernel,model);
       catch exception
           throw(exception)
       end
   end

   
%-------------------------------------------------------------------------
%  Kernel weights   

   lags = 0:T-1;
   w    = zeros(T,1);
   x    = lags/b;           
   
   switch Kernel
       case 'TR'
           TR     = (abs(x) <= 1);
           w(TR)  = 1;
       case 'BT'
           BT     = (abs(x) <= 1);                
           w(BT)  = 1-abs(x(BT));
       case 'PZ'
           PZ1    = (abs(x) >= 0) & (abs(x) <= 1/2);
           PZ2    = (abs(x) >= 1/2) & (abs(x) <= 1);  
           w(PZ1) = 1-6*x(PZ1).^2+6*abs(x(PZ1)).^3;
           w(PZ2) = 2*(1-abs(x(PZ2))).^3;
       case 'TH'
           TH     = (abs(x) <= 1);
           w(TH)  = (1+cos(pi*x(TH)))/2;
       case 'QS'
            argQS = 6*pi*x/5;
            w1    = 3./(argQS.^2);
            w2    = (sin(argQS)./argQS)-cos(argQS);
            w     = w1.*w2;
            w(x == 0) = 1;
       otherwise
           error(message('Invalid kernel'))
   end
   
   PhiHat = w(1)*(V'*V)/T;   
   mSigmahat = PhiHat;       % short-run covariance matrix
   DelHat = w(1)*(V'*V)/T;
   for i = 1:T-p-1
       LagCov = (V(1+i:T-p,:)'*V(1:T-p-i,:))/T;
       PhiHat = PhiHat+w(i+1)*(LagCov+LagCov');   
       DelHat = DelHat+w(i+1)*LagCov';             
   end
   % Recolor prewhitened PhiHat and DelHat:
   if p >0
       C           = eye(N)-sum(reshape([VARcoeffs{:}],N,N,p),3);
       LongRunCov  = (C\PhiHat)/C';     % two-sided long-run covariance
       OneSidedLRV = (C\DelHat)/C';     % one-sided long-run covariance
   elseif p == 0
       LongRunCov  = PhiHat;  
       OneSidedLRV = DelHat;
   end
   
   %-----------------------------------------------------------------------
   % This part is specifically designed for the Simulated Estimator in NLS
   % project
   mDelta1 = mSigmahat - OneSidedLRV';    % \mSigma - \mDelta'
   n = 1;
   BiasCorr_SimulatedEst = mDelta1(n+1:end,1:n);
   %-----------------------------------------------------------------------

%-------------------------------------------------------------------------
%  Data-driven bandwidth based on Andrews (1991)  
function b = getBandwidth(V,Kernel,model)

   [T,N]  = size(V);
   wAlpha = ones(N,1);           % Unit weight for slope coefficients
   
   % Set MLE Options
   options = optimset('fmincon');
   options = optimset(options,'display','off','diagnostics','off',...
                                    'algorithm','sqp','tolCon',1e-7);
                                
   switch model                  % page 835 by Andrews (1991)
       case {'AR1','AR1OLS','AR1MLE'}
           rho = zeros(N,1);     % AR coefficient estimates
           sigmaSq = zeros(N,1); % Residual variances
           
           switch model
               case {'AR1','AR1MLE'}
                   AR1 = arima(1,0,0); % AR(1) model
                   for j = 1:N                  
                       ARfit = estimate(AR1,V(:,j),'print',false,'options',options);
                       rho(j) = ARfit.AR{1};
                       sigmaSq(j) = ARfit.Variance;                
                   end
               case 'AR1OLS'
                   for j = 1:N               
                       v      = V(:,j);
                       vt     = v(2:T);
                       vtLag1 = v(1:T-1);
                       rhoj   = vtLag1\vt;
                       res    = vt-rhoj*vtLag1;
                       rho(j) = rhoj;
                       sigmaSq(j) = res'*res/(T-1);
                   end
           end
           num0 = 4*(rho.^2).*(sigmaSq.^2);
           den  = (sigmaSq.^2)./((1-rho).^4);
       case 'ARMA11'
           rho = zeros(N,1);      % AR coefficient estimates
           psi = zeros(N,1);      % MA coefficient estimates
           sigmaSq = zeros(N,1);  % Residual variances           
           ARMA11 = arima(1,0,1); % ARMA(1,1) model
           
           for j = 1:N          
               v = V(:,j);
               ARMAfit = estimate(ARMA11,v,'print',false,'options',options);
               rho(j) = ARMAfit.AR{1};
               psi(j) = ARMAfit.MA{1};
               sigmaSq(j) = ARMAfit.Variance;
           end
           num0 = 4*((1+rho.*psi).^2).*((rho+psi).^2).*(sigmaSq.^2);
           den = (((1+psi).^4).*(sigmaSq.^2))./((1-rho).^4);
   end
   
   % Compute alphas
   num1 = ((1-rho).^6).*((1+rho).^2);
   num2 = (1-rho).^8;
   denSum = sum(wAlpha.*den);
   alpha1 = sum(wAlpha.*(num0./num1))/denSum;
   alpha2 = sum(wAlpha.*(num0./num2))/denSum;
   
   % Compute optimal bandwidth
   switch Kernel                  % page 834 by Andrews (1991)
       case 'TR'       
           b = (0.6611)*(alpha2*T)^(1/5);        
       case 'BT'       
           b = (1.1447)*(alpha1*T)^(1/3);               
       case 'PZ'       
           b = (2.6614)*(alpha2*T)^(1/5);              
       case 'TH'      
           b = (1.7462)*(alpha2*T)^(1/5);                
       case 'QS'      
           b = (1.3221)*(alpha2*T)^(1/5);       
       otherwise
           error(message('Bandwidth selection: invalid Kernel'))
   end
   
   
%-------------------------------------------------------------------------
%  Create a cell array of NaN's
%  N = m+n; p: number of lags
function NaNp = nancell(N,p)

   NaNp    = cell(1,p);
   NaNp(:) = {NaN(N)};  
   
   
%-------------------------------------------------------------------------
%  Check the validity of 'Bandwidth' input
function OK = CheckBandwidth(b)

   if ~isvector(b)
       error(message('Bandwidth is not a vector.'))     
   elseif isnumeric(b) && ~isscalar(b)
       error(message('Bandwidth value is not a scalar.'))      
   elseif ischar(b) && ~ismember(upper(b),{'AR1','AR1OLS','AR1MLE','ARMA11'})    
       error(message('Bandwidth estimation is not valid.'))
   else
       OK = true;
   end      
         
   
  













