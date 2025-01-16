%########## Multivariate BIAM: Lin and Reuvers (2019, working paper)  ##########%

% This function gives all the covariance estimates that are necessary 
% for applying FMGLS.

% <<<<< List of Inputs >>>>> %
% *1.  U = [vu_1'
%           vu_2'
%            ...
%           vu_T'] is a [T x n] matrix where vu_t is a [n x 1] vector;
%      V is a [T x m] matrix with vv_t a [m x 1] vector.
% *2 'Prewhiten': (1) the default programme does not apply prewhitening
%                 procedure for the process {\vu_t};
%                  (2) if the user would like to apply this procedure, 
%                 specify ['Prewhiten','VAR']. Then it can prewhiten the 
%                 {\vu_t} process by a VAR(1) process.
% *3 'SARNorm': the matrix norm that is used for selecting optimal banding
%               parameter, please specify one of {1,2,'inf','fro'}.


% <<<<< List of Outputs >>>>> %
% Three quantities are reporte; check "estimation_procedure.pdf".
%
% Q1: BIAM \widehat{\mOmega}_{u,nT}^{-1};
% Q2: \widehat{\mOmega}_{vv}^{-1} \widehat{\mOmega}_{vu} \widehat{\mOmega}_{uu}^{-1};
% Q3: \widehat{\mDelta}_{vv}
% Q4: \widehat{\mOmega}_{uGv}, which is used for KPSS and Wald tests
% Q5: \mSigma
% Q6: \widehat{\mOmega}_{xi}, the long-run vairance of [\vu_t,\vvt]'
% Q7: {\widehat{\mSigma}-\widehat{\mDelta}'}_{vu}, for bias correction in NLS
% Q8: \widehat{\mOmega}_{uu}^{-1} \widehat{\mOmega}_{uGv} \widehat{\mOmega}_{uu}^{-1} for Wald tests
% Q9: \widehat{\mOmega}_{uu}^{-1}

% <<<<< Examples of Executing this Function >>>>> %
%  case 1, no prewhitening, default matrix norm for banding parameter
%          selection to be L_1 norm, implement
%    [Q1,Q2,Q3,Q4,Q5,Q6,Q7] = LRbiam(U,V);
% 
%  For feasible FMGLS, only the first three quantities are used, run
%    [Q1,Q2,Q3] = LRbiam(U,V);
%  is sufficient.
% 
%  case 2, with prewhitening
%    [Q1,Q2,Q3,Q4] = LRbiam(U,V,'Prewhiten','VAR');
% 
%  case 3, Normally, the band selection is not sensitive to matrix norm. If
%          the user would like to change the matrix norm for band selection, 
%          say Frobenius norm, try
%    [Q1,Q2,Q3,Q4] = LRbiam(U,V,'SARNorm','fro');
          


function [Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9] = LRbiam(U,V,varargin)

     [T,n] = size(U);      % T: sample size; n: dimension of \vu_t
     [~,m] = size(V);      % m: dimension of \vv_t
     Xi    = [U V];
     N     = n+m;
    
%<<<<< Parse inputs and set defaults >>>>>%

     defaultPrewhiten  = 'NONE';       
     expectedPrewhiten = {'NONE','VAR'}; 
     defaultToEstBIAM  = 'ON';       
     expectedToEstBIAM = {'ON','OFF'};      
     defaultSARNorm    = 1;       
     
     p = inputParser;
     addParameter(p,'Prewhiten',defaultPrewhiten,@(x) any(validatestring(x,expectedPrewhiten))); 
     addParameter(p,'ToEstBIAM',defaultToEstBIAM ,@(x) any(validatestring(x,expectedToEstBIAM )));      
     addParameter(p,'SARNorm',defaultSARNorm,@Checksarnorm);
     parse(p,varargin{:});
     
     prewhiten = p.Results.Prewhiten;
     ToEstBIAM = p.Results.ToEstBIAM;     
     bnorm     = lower(p.Results.SARNorm);
     
     
%-------------------------------------------------------------------------
%  Prewhitening procedure based on VAR(1) & multivariate BIAM for {\vu_t}

     switch ToEstBIAM
         case 'ON'
             switch prewhiten 
                 case 'NONE'                          % non-prewhitened    
                     Utran = U'; 
                     Q1    = bimam(T,n,Utran,bnorm);  % multivariate BIAM estimate of {\vu_t}
                 case 'VAR'                           % prewhitened by VAR(1)
                     U1    = U';
                     A1    = (U1(:,2:T)*U1(:,1:T-1)')/(U1(:,1:T-1)*U1(:,1:T-1)'); % VAR(1) OLS Regression
                     D1    = kron(eye(T),eye(n))-kron(diag(ones(T-1,1),-1),A1);
                     Utran = reshape(D1*U1(:),[n,T]);
                     Q1w   = bimam(T,n,Utran,bnorm);      % multivariate BIAM estimate of {\vu_t} before recoloring
                     Q1    = D1'*Q1w *D1;                 % recoloring
             end
         case 'OFF'
             Q1 = [];
     end


     
%-------------------------------------------------------------------------     
%  Estimating Quantities Q2 - Q9
 
     [mOmega,mSig,Dvv,mDelta] = multiAM(T,n,m,N,Xi,bnorm);
        
     % \widehat{\mOmega}_{vv}^{-1} \widehat{\mOmega}_{vu} \widehat{\mOmega}_{uu}^{-1}
     mOmega_uu = mOmega(1:n,1:n);
     mOmega_uv = mOmega(1:n,end-m+1:end);
     mOmega_vv = mOmega(end-m+1:end,end-m+1:end);
     Q2        = mOmega_vv\mOmega_uv'/mOmega_uu;     
     
     % \mDelta_vv: one-sided LRV
     Q3 = Dvv;     
     
     % \widehat{\mOmega}_{uGv}
     Q4 = mOmega_uu-(mOmega_uv/mOmega_vv)*mOmega_uv';
     
     % \mSigma
     Q5 = mSig;
     
     % \widehat{\mOmega}_{xi}
     Q6 = mOmega;
     
     % \widehat{\mOmega}_{zeta}
     mDelta1 = mSig-mDelta';        % \mSigma - \mDelta'
     Q7 = mDelta1(n+1:end,1:n);
     
     % \widehat{\mOmega}_{uu}^{-1} \widehat{\mOmega}_{uGv} \widehat{\mOmega}_{uu}^{-1}
     Q8 = mOmega_uu\Q4/mOmega_uu;
     
     % \widehat{\mOmega}_{uu}^{-1}
     Q9 = inv(mOmega_uu);
     
%-------------------------------------------------------------------------
%  Special design: function that provides estimates for Q2 - Q5
%  mOmega: Long-run covariance matrix
%  mSig:   \mSigma, the variance of [\veta_t',\vepsi_t']'
%  Dvv:    one-sided long run covariance of {\vv_t}
%  mDelta: one-sided long run covariance of {\vxi_t}

function [mOmega,mSig,Dvv,mDelta] = multiAM(T,n,m,N,Xi,bnorm)

    Xitran = Xi';
    qT     = sar(T,N,Xitran,bnorm);   % choose the banding parameter by SAR
    
    % banded autocovariance matrix of {\vxi_t};  not its inverse
    Mq = kron(eye(T),eye(N));   
    Sq = kron(eye(T),zeros(N,N)); 
    Sq(1:N,1:N) = Xitran*Xitran'/T; 
    for j = 1:qT
        yy     = Xitran(:,j+1:end);
        xxlagj = lagmatrix(Xitran',1:1:j)';         % taking lags
        xx     = xxlagj(:,j+1:end);                 % lag regressors
        Aj     = yy*xx'/(xx*xx');
        
        ordercell = cellfun(@(x,y)x:y,num2cell([flip(1:j-1) 0]*N+1),...
                             num2cell(flip(1:j)*N),'UniformOutput',false);
        reorder   = [ordercell{:}];
        
        Mq((j*N+1):(j+1)*N,1:(j*N)) = -Aj(:,reorder);
        ej = yy - Aj*xx;
        Sq((j*N+1):(j+1)*N,(j*N+1):(j+1)*N) = (ej*ej')/(size(ej,2));
    end   
    Aq  = Aj(:,reorder); 
    SNq = (ej*ej')/(size(ej,2));
    for j = qT+1:T-1
       Mq((j*N+1):(j+1)*N,(j-qT)*N+1:(j*N)) = -Aq;  
       Sq((j*N+1):(j+1)*N,(j*N+1):(j+1)*N)  = SNq;
    end 
    OmegaXi = Mq\Sq/Mq';   
    
    % mSig
    mSig = SNq;
        
    % Long-run covariance matrix of \vxi_t
    H = eye(N);                   % to collect the autoregressive coefficient matrices
    minusAj = -Aj;
    for j = 1:qT
       jk = (j-1)*N+1:j*N;
       H = H + minusAj(:,jk);
    end
    A1     = H(1:n,1:n);
    invD1  = H(end-m+1:end,end-m+1:end);
    H1     = blkdiag(A1,invD1);    % H1 = [ A1   0
                                   %        0  invD1]
    mOmega = H1\mSig/(H1');        %  Long-run covariance matrix of \vzeta_t
    
    % Dvv
    SelV   = [zeros(m,n) eye(m)];         % selection matrix to choose \vv_t from \vxi_t
    SS     = kron(eye(T),SelV);
    OmegaV = SS*OmegaXi*SS';              % banded autocovariance matrix of {\vv_t};
    i      = (T-1)*m+1:T*m;
    Dvv    = zeros(m,m);
    dr     = 0.01;
    rT     = min(ceil(0.5*T/qT^(3+dr)),T);
    for j = (T-rT+1):T
       jk  = (j-1)*m+1:j*m;
       Dvv = Dvv + OmegaV(i,jk)';
    end
    
    % mDelta_xi
    i      = (T-1)*(m+n)+1:T*(m+n);
    mDelta = zeros(m+n,m+n);
    for j = (T-rT+1):T
       jk  = (j-1)*(m+n)+1:j*(m+n);
       mDelta = mDelta + OmegaXi(i,jk)';
    end
    
    
%-------------------------------------------------------------------------
%  the multivariate banded inverse autocovariance matrix   
function Mult_BIAM = bimam(T,n,Utran,bnorm)
  
    qT = sar(T,n,Utran,bnorm);   % choose the banding parameter by SAR

    Mq    = kron(eye(T),eye(n));   
    Sqinv = kron(eye(T),zeros(n,n));
    Sqinv(1:n,1:n) = inv(Utran*Utran'/T);
    for j = 1:qT
        yy     = Utran(:,j+1:end);
        xxlagj = lagmatrix(Utran',1:1:j)';         % taking lags
        xx     = xxlagj(:,j+1:end);                % lag regressors
        Aj     = yy*xx'/(xx*xx');
        
        ordercell = cellfun(@(x,y)x:y,num2cell([flip(1:j-1) 0]*n+1),...
                             num2cell(flip(1:j)*n),'UniformOutput',false);
        reorder   = [ordercell{:}];
        
        Mq((j*n+1):(j+1)*n,1:(j*n)) = -Aj(:,reorder);
        ej = yy - Aj*xx;
        Sqinv((j*n+1):(j+1)*n,(j*n+1):(j+1)*n) = inv((ej*ej')/(size(ej,2)));
    end
    Aq     = Aj(:,reorder);
    Snqinv = inv((ej*ej')/(size(ej,2)));
    for j = qT+1:T-1
       Mq((j*n+1):(j+1)*n,(j-qT)*n+1:(j*n))   = -Aq;  
       Sqinv((j*n+1):(j+1)*n,(j*n+1):(j+1)*n) = Snqinv;
    end
    Mult_BIAM = Mq'*Sqinv*Mq;
    
    
%-------------------------------------------------------------------------
%  selection of the banding length q 
%  Wtran: any matrix that has similar construction as Utran above
function [qT,RiskVal_min] = sar(T,n,Wtran,bnorm)

    warning('off','all')

    ql = 1;                                      
    qu = ceil(2*T^(0.25));                     % \bar{q}
    w0 = floor(T/5);                           % length of subseries
    N0 = floor(T/w0);                          % number of subseries   
    
    % the inverse of sample multivariate autocovariance matrix Pi
    Pi_inv = smamInv(T,n,Wtran,qu);
    
    % selecting the banding parameter qT
    Usplit  = reshape(Wtran(:,1:w0*N0),[n*w0,N0]);  % splict residuals into N0
                                                    % subsequences (every column in the matrix)
    RiskVal = zeros(qu-ql,1);                       % to store the risk values
    count   = 0;
    for b = ql:qu-1
        count = count+1;
        rv    = zeros(N0,1);
        for n0 = 1:N0
            Un0   = reshape(Usplit(:,n0),[n,w0]);
            Mq    = kron(eye(qu),eye(n)); 
            Sqinv = kron(eye(qu),zeros(n,n));
            Sqinv(1:n,1:n) = inv(Un0*Un0'/T);
            for j = 1:b
                yy     = Un0(:,j+1:end);
                xxlagj = lagmatrix(Un0',1:1:j)';         % taking lags
                xx     = xxlagj(:,j+1:end);             % lag regressors
                Anj    = yy*xx'/(xx*xx');
                
                ordercell = cellfun(@(x,y)x:y,num2cell([flip(1:j-1) 0]*n+1),...
                             num2cell(flip(1:j)*n),'UniformOutput',false);
                reorder   = [ordercell{:}];
                
                Mq((j*n+1):(j+1)*n,1:(j*n)) = -Anj(:,reorder);
                ej = yy - Anj*xx;
                Sqinv((j*n+1):(j+1)*n,(j*n+1):(j+1)*n) = inv((ej*ej')/(size(ej,2)));
            end
            Anq    = Anj(:,reorder);
            Snqinv = inv((ej*ej')/(size(ej,2)));
            for j = b+1:qu-1
                Mq((j*n+1):(j+1)*n,(j-b)*n+1:(j*n))    = -Anq;  
                Sqinv((j*n+1):(j+1)*n,(j*n+1):(j+1)*n) = Snqinv;
            end
            bimam  = Mq'*Sqinv*Mq;
            D      = bimam - Pi_inv;
            rv(n0) = norm(D,bnorm);
        end
        RiskVal(count) = N0^(-1)*sum(rv);  
    end
    [RiskVal_min,m_min] = min(RiskVal);
    qT = m_min+ql-1;    

    
%-------------------------------------------------------------------------
%  the inverse of qu*qu sample multivariate autocovariance matrix Pi
function Pi_inv = smamInv(T,n,Wtran,qu)

    Pi = zeros(n*qu,n*qu);
    for t = qu:T-1
        Et = reshape(Wtran(:,t:-1:t-qu+1),[n*qu,1]);
        Pi = Pi + Et*Et';
    end
    Pi = (T-qu)^(-1)*Pi;
    Pi_inv = inv(Pi); 


%-------------------------------------------------------------------------
%  Check the validity of 'SARNorm' input
function OK = Checksarnorm(mnorm)

   if isnumeric(mnorm) && mnorm ~= 1 && mnorm ~=2
       error('Matrix norm must take 1, 2, ''inf'' or ''fro''.')  
   elseif ischar(mnorm) && ~ismember(lower(mnorm),{'inf','fro'})
       error('Matrix norm must take 1, 2, ''inf'' or ''fro''.')  
   else 
       OK = true;
   end
   
  
