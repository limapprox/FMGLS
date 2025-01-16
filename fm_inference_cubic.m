%###################### Fully Modified Estimation #######################%

function [vbetaFMOLS,vbetaFMSUR,vbetaFMGLS,K_ct,WaldstatLinear,WaldstatQuadra,WaldstatCubic,Table] = fm_inference_cubic(n,T,vy,mZ,vx,Vhat,varargin)

defaultOrganizeOutputs = 'NO';                  
expectedOrganizeOutputs  = {'NO','YES'}; 
defaultSelectedCountries  = 'NONE';                  

p = inputParser;
addOptional(p,'OrganizeOutputs',defaultOrganizeOutputs,@(x) any(validatestring(x,expectedOrganizeOutputs))); 
addOptional(p,'SelectedCountries',defaultSelectedCountries);
parse(p,varargin{:});

OrganizeOutputs = upper(p.Results.OrganizeOutputs);
countries_selected = p.Results.SelectedCountries;

VhatTran = Vhat';
vv   = VhatTran(:);
Npar = 4*n;

%===================== FIRST-STEP HAC/BIAM QUANTITIES ====================%
% 1. HAC
vbeta_ols = (mZ'*mZ)\(mZ'*vy);
vuhat = vy - mZ*vbeta_ols; 
Uhat  = reshape(vuhat,[n,T])';
[mOmegahat_uvHAC,mOmegahat_vvHAC,mOmegahat_uGvHAC,mDeltahat_vuplusHAC,PreMathatHAC] = hac_quantities(n,T,Uhat,Vhat);

% 2. BIAM
% [invCovhat,Q2,mDeltahat_vvBIAM,mOmegahat_uGvBIAM,mSigmahatBIAM,mOmegahatBIAM,~,Q8,Q9] = LRbiam(Uhat,Vhat);
[invCovhat,Q2,mDeltahat_vvBIAM,~,mSigmahatBIAM,~,~,Q8,Q9] = LRbiam(Uhat,Vhat);
mSigmahat_epsetaBIAM = mSigmahatBIAM(n+1:end,1:n);
mSigmahat_etaetaBIAM = mSigmahatBIAM(1:n,1:n);
%=========================================================================%

%============================= FM ESTIMATIONS ============================%
% 1. FM-SOLS
vyplus_fOLS = vy - kron(eye(T),mOmegahat_uvHAC/mOmegahat_vvHAC)*vv;
vA_fOLS = zeros(Npar,1);
for i = 1:n  
    vA_fOLS(4*(i-1)+1:4*i) = mDeltahat_vuplusHAC(i,i)*[zeros(1,1); T; 2*sum(vx(i,:)); 3*sum(vx(i,:).^2)];       
end
vbetaFMOLS = (mZ'*mZ)\(mZ'*vyplus_fOLS-vA_fOLS);

% 2. FM-SUR
vyplus_fSUR = vyplus_fOLS;
vA_fSURcoef = mDeltahat_vuplusHAC/mOmegahat_uGvHAC;   % the coefficient matrix in Equation (12) of WGH(2019) 
vA_fSUR     = zeros(Npar,1);                          % \tilde{\vA}^* vector in Equation (12) 
for i = 1:n
    vA_fSUR(4*(i-1)+1:4*i) = vA_fSURcoef(i,i)*[zeros(1,1); T; 2*sum(vx(i,:)); 3*sum(vx(i,:).^2)];    
end
vbetaFMSUR = (mZ'*PreMathatHAC*mZ)\(mZ'*PreMathatHAC*vyplus_fSUR-vA_fSUR);

% 3. FM-GLS
vB_fGLScoef = mSigmahat_epsetaBIAM/mSigmahat_etaetaBIAM-mDeltahat_vvBIAM*Q2;
vB_fGLS = zeros(Npar,1);                        
for i = 1:n
    vB_fGLS(4*(i-1)+1:4*i) = vB_fGLScoef(i,i)*[zeros(1,1); T; 2*sum(vx(i,:)); 3*sum(vx(i,:).^2)]; 
end
PreMathatGLS = kron(eye(T),Q2');
vbetaFMGLS = (mZ'*invCovhat*mZ)\(mZ'*invCovhat*vy-mZ'*PreMathatGLS*vv-vB_fGLS);
%=========================================================================%

%============================== Nonlinear CT =============================%
% K_ct = zeros(4,4);     % SOLS | SUR | FGLS | BIAM
K_ct = zeros(4,3);     % SOLS | SUR | BIAM

% 1. KPSS-SOLS
vuhatplus_SOLS = vyplus_fOLS - mZ*vbetaFMOLS;
[Kstat_SOLS,rej_rule_SOLS,bopt_SOLS,Mopt_SOLS] = multikpss_bonferroni(T,n,vuhatplus_SOLS,PreMathatHAC);
K_ct(:,1) = [Kstat_SOLS,rej_rule_SOLS,bopt_SOLS,Mopt_SOLS]';

% 2. KPSS-SUR
vuhatplus_SUR = vyplus_fSUR - mZ*vbetaFMSUR;
[Kstat_SUR,rej_rule_SUR,bopt_SUR,Mopt_SUR] = multikpss_bonferroni(T,n,vuhatplus_SUR,PreMathatHAC);
K_ct(:,2) = [Kstat_SUR,rej_rule_SUR,bopt_SUR,Mopt_SUR];

% 3. KPSS-FGLS
% mOmegahat_uvBIAM = mOmegahatBIAM(1:n,n+1:end); 
% mOmegahat_vvBIAM = mOmegahatBIAM(n+1:end,n+1:end);   
% vuhatplus_fGLS = vy - kron(eye(T),mOmegahat_uvBIAM/mOmegahat_vvBIAM)*vv - mZ*vbetaFMGLS;       
% [Kstat_fGLS,rej_rule_fGLS,bopt_fGLS,Mopt_fGLS] = multikpss_bonferroni(T,n,vuhatplus_fGLS,kron(eye(T),inv(mOmegahat_uGvBIAM))); 
% K_ct(:,3) = [Kstat_fGLS,rej_rule_fGLS,bopt_fGLS,Mopt_fGLS]';

% 3. KPSS-BIAM
vuhat_fGLS = vy - mZ*vbetaFMGLS; 
[Kstat_BIAM,rej_rule_BIAM,bopt_BIAM,Mopt_BIAM] = multikpss_bonferroni(T,n,vuhat_fGLS,invCovhat); 
K_ct(:,3) = [Kstat_BIAM,rej_rule_BIAM,bopt_BIAM,Mopt_BIAM]';
%=========================================================================%

%=============================== Wald Tests ==============================%
mRlinear = kron(eye(n),[0 1 0 0]);    % select the linear coefficient
mRquadra = kron(eye(n),[0 0 1 0]);    % select the quadratic coefficient
mRcubic = kron(eye(n),[0 0 0 1]);     % select the cubic coefficient
WaldstatLinear = zeros(n,3);    % SOLS | SUR | FGLS 
WaldstatQuadra = zeros(n,3);    % SOLS | SUR | FGLS
WaldstatCubic = zeros(n,3);    % SOLS | SUR | FGLS

% 1. SOLS
mPhi_FMOLSlinear = mRlinear*((mZ'*mZ)\(mZ'*kron(eye(T),mOmegahat_uGvHAC)*mZ)/(mZ'*mZ))*mRlinear'; 
mPhi_FMOLSquadra = mRquadra*((mZ'*mZ)\(mZ'*kron(eye(T),mOmegahat_uGvHAC)*mZ)/(mZ'*mZ))*mRquadra'; 
mPhi_FMOLScubic = mRcubic*((mZ'*mZ)\(mZ'*kron(eye(T),mOmegahat_uGvHAC)*mZ)/(mZ'*mZ))*mRcubic'; 
for i = 1:n
    
    mR1 = mRlinear(i,:);
    mR2 = mRquadra(i,:);
    mR3 = mRcubic(i,:);
    
    WaldstatLinear(i,1) = (mR1*vbetaFMOLS)^2/mPhi_FMOLSlinear(i,i); % Single-equation Wald test
    WaldstatQuadra(i,1) = (mR2*vbetaFMOLS)^2/mPhi_FMOLSquadra(i,i); % Single-equation Wald test
    WaldstatCubic(i,1) = (mR3*vbetaFMOLS)^2/mPhi_FMOLScubic(i,i); % Single-equation Wald test
  
end

% 2. SUR
mPhi_FMSURlinear = (mRlinear/(mZ'*PreMathatHAC*mZ))*mRlinear';
mPhi_FMSURquadra = (mRquadra/(mZ'*PreMathatHAC*mZ))*mRquadra';
mPhi_FMSURcubic = (mRcubic/(mZ'*PreMathatHAC*mZ))*mRcubic';
for i = 1:n
    
    mR1 = mRlinear(i,:);
    mR2 = mRquadra(i,:);
    mR3 = mRcubic(i,:);
    
    WaldstatLinear(i,2) = (mR1*vbetaFMSUR)^2/mPhi_FMSURlinear(i,i); % Single-equation Wald test
    WaldstatQuadra(i,2) = (mR2*vbetaFMSUR)^2/mPhi_FMSURquadra(i,i); % Single-equation Wald test
    WaldstatCubic(i,2) = (mR3*vbetaFMSUR)^2/mPhi_FMSURcubic(i,i); % Single-equation Wald test
    
end

% 3. FGLS
mPhi_FMGLSlinear = mRlinear*((mZ'*kron(eye(T),Q9)*mZ)\(mZ'*kron(eye(T),Q8)*mZ)/(mZ'*kron(eye(T),Q9)*mZ))*mRlinear';
mPhi_FMGLSquadra = mRquadra*((mZ'*kron(eye(T),Q9)*mZ)\(mZ'*kron(eye(T),Q8)*mZ)/(mZ'*kron(eye(T),Q9)*mZ))*mRquadra';
mPhi_FMGLScubic = mRcubic*((mZ'*kron(eye(T),Q9)*mZ)\(mZ'*kron(eye(T),Q8)*mZ)/(mZ'*kron(eye(T),Q9)*mZ))*mRcubic';
for i = 1:n
    
    mR1 = mRlinear(i,:);
    mR2 = mRquadra(i,:);
    mR3 = mRcubic(i,:);
    
    WaldstatLinear(i,3) = (mR1*vbetaFMGLS)^2/mPhi_FMGLSlinear(i,i); % Single-equation Wald test
    WaldstatQuadra(i,3) = (mR2*vbetaFMGLS)^2/mPhi_FMGLSquadra(i,i); % Single-equation Wald test
    WaldstatCubic(i,3) = (mR3*vbetaFMGLS)^2/mPhi_FMGLScubic(i,i); % Single-equation Wald test
    
end
%=========================================================================%

%============================ Organize Outputs ===========================%

switch OrganizeOutputs
case 'NO'
    Table = 'Empty';
case 'YES'
    OutputTable = zeros(3*n,3);   %         |            | SOLS | SUR | FGLS
                                  %         | beta_{i,1} |
                                  % country | beta_{i,2} |
                                  %         | beta_{i,3} |
                                  %---------------------------------------------
    for i = 1:n
        
        mR1 = mRlinear(i,:);
        mR2 = mRquadra(i,:);
        mR3 = mRcubic(i,:);
        
        beta1FMOLS = mR1*vbetaFMOLS;
        beta2FMOLS = mR2*vbetaFMOLS;
        beta3FMOLS = mR3*vbetaFMOLS;
        OutputTable(3*(i-1)+(1:3),1) = [round(beta1FMOLS,3); round(beta2FMOLS,3); round(1E5*beta3FMOLS,3)];
        
        beta1FMSUR = mR1*vbetaFMSUR;
        beta2FMSUR = mR2*vbetaFMSUR;
        beta3FMSUR = mR3*vbetaFMSUR;
        OutputTable(3*(i-1)+(1:3),2) = [round(beta1FMSUR,3); round(beta2FMSUR,3); round(1E5*beta3FMSUR,3)];
        
        beta1FMGLS = mR1*vbetaFMGLS;
        beta2FMGLS = mR2*vbetaFMGLS;
        beta3FMGLS = mR3*vbetaFMGLS;
        OutputTable(3*(i-1)+(1:3),3) = [round(beta1FMGLS,3); round(beta2FMGLS,3); round(1E5*beta3FMGLS,3)];
    end
    Table = table(reshape(repmat(countries_selected,[3,1]),[3*n,1]),OutputTable);
end

%=========================================================================%



end

function [mOmegahat_uvHAC,mOmegahat_vvHAC,mOmegahat_uGvHAC,mDeltahat_vuplusHAC,PreMathatHAC] = hac_quantities(n,T,Uhat,Vhat)

[mOmegahatHAC,mDeltahatHAC] = AndHAC([Uhat Vhat],'Kernel','BT');
mOmegahat_uuHAC = mOmegahatHAC(1:n,1:n);
mOmegahat_vvHAC = mOmegahatHAC(n+1:end,n+1:end);
mOmegahat_vuHAC = mOmegahatHAC(n+1:end,1:n);
mOmegahat_uvHAC = mOmegahatHAC(1:n,n+1:end);
mDeltahat_vuHAC = mDeltahatHAC(n+1:end,1:n);
mDeltahat_vvHAC = mDeltahatHAC(n+1:end,n+1:end);
mDeltahat_vuplusHAC = mDeltahat_vuHAC-mDeltahat_vvHAC*(mOmegahat_vvHAC\mOmegahat_vuHAC);
mOmegahat_uGvHAC = mOmegahat_uuHAC-mOmegahat_uvHAC*(mOmegahat_vvHAC\mOmegahat_vuHAC);
PreMathatHAC = kron(eye(T),inv(mOmegahat_uGvHAC));

end






