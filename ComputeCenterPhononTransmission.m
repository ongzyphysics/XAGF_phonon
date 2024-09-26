function PhononData = ComputeCenterPhononTransmission(w_in,ww_in,LeftLead,RightLead,Center,Left,Right)
% NOTE: This function calculates the phonon transmission for frequency range wvec
% as given. LeftLead and RightLead contain the surface Greens functions and
% the modes in the leads. Left, Center and Right contain the force constant
% matrices.
% NOTE: This function requires "GreensFunctionSolver.m"

  Lyr = Center.Lyr;   % data structs containing HL, HC, HR submatrices describing Hamiltonian of scattering region
  nmax = length(Lyr); % number of layers in scattering region

  HCL = Center.HCL; % coupling matrix from scattering region to left lead
  HLC = HCL';       % hermitian conjugate of above 
  HCR = Center.HCR; % coupling matrix from scattering region to right lead 
  HRC = HCR';       % hermitian conjugate of above 

  a_L = Left.a_long;  % longitudinal lattice constant (left lead)
  a_R = Right.a_long; % longitudinal lattice constant (right lead)
  MatHC_LeftMost  = Lyr(1).MatHC; 
  MatHC_RightMost = Lyr(nmax).MatHC;

  clearvars Center Left Right

  % === Pre-allocate PhononData ===
  PhononData.Xi_mode     = 0; 
          % total transmission by summing over all channel transmission in right lead
  PhononData.AntiXi_mode = 0; 
          % total reflection summing over all channel reflection in left lead 
  PhononData.Xi_negf     = 0; 
          % total transmission from Caroli's formula

  % ===== (t and r submatrices for left-to-right and right-to-left fluxes) =====
  %{
  PhononData.t_RL  = zeros(size(RightLead(1).U_plus,1),size(LeftLead(1).U_minus_adv,1));
          % t-matrix for transition amplitudes from left-lead to right-lead channels 
  PhononData.t_LR  = zeros(size(LeftLead(1).U_minus,1),size(RightLead(1).U_plus_adv,1));
          % t-matrix for transition amplitudes from right-lead to left-lead channels 
  PhononData.r_LL  = zeros(size(LeftLead(1).U_minus,1),size(LeftLead(1).U_minus_adv,1));
          % r-matrix for transition amplitudes of left-lead backscattered/reflected channels
  PhononData.r_RR  = zeros(size(RightLead(1).U_plus,1),size(RightLead(1).U_plus_adv,1));
          % r-matrix for transition amplitudes of right-lead backscattered/reflected channels
  %}
  
  % ===== (transmission, absorption and reflection coefficients) =====
  PhononData.xi_inc_L = diag(zeros(size(LeftLead(1).U_minus_adv)));
          % matrix of transmission coefficients of left-lead channels (left-to-right flux)
  PhononData.xi_abs_L = diag(zeros(size(RightLead(1).U_plus)));
          % matrix of absorption coefficients of right-lead channels (left-to-right flux)
  PhononData.xi_ref_L = diag(zeros(size(LeftLead(1).U_minus)));
          % matrix of reflection coefficients of left-lead channels (left-to-right flux)

  PhononData.xi_inc_R = diag(zeros(size(RightLead(1).U_plus_adv)));
          % matrix of transmission coefficients of right-lead channels (right-to-left flux)
  PhononData.xi_abs_R = diag(zeros(size(LeftLead(1).U_minus)));
          % matrix of absorption coefficients of left-lead channels (right-to-left flux)
  PhononData.xi_ref_R = diag(zeros(size(RightLead(1).U_plus)));
          % matrix of reflection coefficients of right-lead channels (right-to-left flux)

  %{
  PhononData.t    = zeros(size(RightLead(1).U_plus,1),size(LeftLead(1).U_minus_adv,1));
          % t-matrix for transition amplitudes from left-lead to right-lead channels 
  PhononData.tt_R = diag(zeros(size(RightLead(1).U_plus)));
          % matrix of transmission coefficients of right-lead channels (transmitted)
  PhononData.tt_L = diag(zeros(size(LeftLead(1).U_minus_adv)));
          % matrix of transmission coefficients of left-lead channels (incident)
  PhononData.r    = zeros(size(LeftLead(1).U_minus,1),size(LeftLead(1).U_minus_adv,1));
          % matrix for transition amplitudes of left-lead backscattered/reflected channels
  PhononData.rr_L = zeros(size(LeftLead(1).U_minus,1),size(LeftLead(1).U_minus,1));
          % matrix of reflection coefficients of lead-lead channels (reflected)
  %}

  % === Compute transmission and reflection ===

  % === EFFECTIVE HAMILTONIAN AND GREEN'S FUNCTION OF DECOUPLED SCATTERING REGION === 
  SurfGR = RightLead.MatSurfGR; % right-lead surface GF
  SurfGL = LeftLead.MatSurfGL;  % left-lead surface GF
  % Sigma_R = Center.HCR*RightLead.MatSurfGR*(Center.HCR)'; % self-energy matrix corresponding to interaction with right lead 
  % Sigma_L = Center.HCL*LeftLead.MatSurfGL*(Center.HCL)'; % self-energy matrix corresponding to interaction with left lead  

  ww = ww_in; % frequency-square point
  w  = w_in;  % frequency 

  if gt(nmax,1) % if more than one layer in the scattering region
      Lyr(1).MatHC    = MatHC_LeftMost  + HCL*SurfGL*HLC; 
      % Lyr(1).MatHC    = Lyr(1).MatHC  + Sigma_L; 
      % add contribution from left-lead surface GF to get effective Hamiltonian of decoupled scattering region
      Lyr(nmax).MatHC = MatHC_RightMost + HCR*SurfGR*HRC;
      % Lyr(nmax).MatHC = Lyr(nmax).MatHC + Sigma_R;
      % add contribution from right-lead surface GF to get effective Hamiltonian of decoupled scattering region
  else % if only one layer, won't really work since the minimum no. of layers is now three, i.e. nmax >= 3
      Lyr(1).MatHC = MatHC_LeftMost  + HCL*SurfGL*HLC + HCR*SurfGR*HRC;
      % Lyr(1).MatHC = Lyr(1).MatHC  + Sigma_L + Sigma_R;
  end

  [EffG_LR, EffG_RR, EffG_RL, EffG_LL, P] = GreensFunctionSolver(ww,Lyr); 
  % get submatrices in retarded GF of decoupled scattering region   

  clearvars Lyr;
    
  % === TRANSITION AMPLITUDE CALCULATIONS FOR FORWARD AND BACKWARD SCATTERING ===
  
  InvU_plus_R = pinv(RightLead.U_plus);                    % inverse of eigenmode matrix for out-going channels in right lead
  SqrtV_plus_R = sqrt(diag(RightLead.VecV_plus));          % sqrt of eigenvelocity matrix for out-going channels in right lead
  InvU_minus_L = pinv(LeftLead.U_minus);                   % inverse of eigenmode matrix for out-going channels in left lead
  SqrtV_minus_L = sqrt(diag(LeftLead.VecV_minus));         % sqrt of eigenvelocity matrix for out-going channels in left lead

  InvTranU_minus_L_adv = pinv(LeftLead.U_minus_adv');      % inverse-transpose of eigenmode matrix for in-coming channels in left lead
  SqrtV_minus_L_adv = sqrt(diag(LeftLead.VecV_minus_adv)); % sqrt of eigenvelocity matrix for in-coming channels in left lead
  InvTranU_plus_R_adv = pinv(RightLead.U_plus_adv');       % inverse-transpose of eigenmode matrix for in-coming channels in right lead
  SqrtV_plus_R_adv = sqrt(diag(RightLead.VecV_plus_adv));  % sqrt of eigenvelocity matrix for in-coming channels in right lead
  
  EffG0_LL = LeftLead.MatBulkG;
  EffG0_RR = RightLead.MatBulkG;
  
  clearvars LeftLead RightLead
  
  % ===== LEFT-TO-RIGHT ===
  % ===== (transmission matrix) ===
  %{
  U_plus_R = RightLead.U_plus;                   % eigenmode matrix for out-going channels in right lead
  V_plus_R = diag(RightLead.VecV_plus);          % eigenvelocity matrix for out-going channels in right lead
  U_minus_L_adv = LeftLead.U_minus_adv;          % eigenmode matrix for in-coming channels in left lead
  V_minus_L_adv = diag(LeftLead.VecV_minus_adv); % eigenvelocity matrix for in-coming channels in left lead
  t_RL = 2*1i*w*sqrt(V_plus_R)*pinv(U_plus_R)*EffG_RL* ... 
      pinv(U_minus_L_adv')*sqrt(V_minus_L_adv)/sqrt(a_L*a_R);
  %}
  t_RL = 2*1i*w*SqrtV_plus_R*InvU_plus_R*EffG_RL* ... 
      InvTranU_minus_L_adv*SqrtV_minus_L_adv/sqrt(a_L*a_R);
  
  % ===== (reflection matrix) ===
  %{
  U_minus_L = LeftLead.U_minus;          % eigenmode matrix for reflected, out-going channels in left lead
  V_minus_L = diag(LeftLead.VecV_minus); % eigenvelocity matrix for reflected out-going channels in left lead   
  EdgeG  = EffG_LL;
  EdgeG0 = LeftLead.MatBulkG;
  r_LL = 2*1i*w*sqrt(V_minus_L)*pinv(U_minus_L)*(EdgeG-EdgeG0)* ...
      pinv(U_minus_L_adv')*sqrt(V_minus_L_adv)/sqrt(a_L*a_L);
  %}
  r_LL = 2*1i*w*SqrtV_minus_L*InvU_minus_L*(EffG_LL-EffG0_LL)* ...
      InvTranU_minus_L_adv*SqrtV_minus_L_adv/sqrt(a_L*a_L);
  
  % clearvars U_plus_R V_plus_R U_minus_L_adv V_minus_L_adv U_minus_L V_minus_L EdgeG EdgeG0;

  % PhononData.t_RL = t_RL;
  PhononData.t_RL = sparse( t_RL );
  PhononData.Xi_mode = sum(diag(t_RL*t_RL'));
  PhononData.xi_abs_L = diag(t_RL*t_RL');
  PhononData.xi_inc_L = diag(t_RL'*t_RL);

  % PhononData.r_LL = r_LL;
  PhononData.r_LL = sparse( r_LL );
  PhononData.AntiXi_mode = sum(diag(r_LL*r_LL'));
  PhononData.xi_ref_L    = diag(r_LL*r_LL');
  
  clearvars t_RL r_LL;

  %{
  PhononData.t_L = t_L;
  PhononData.Xi_mode = sum(diag(t_L*t_L'));
  PhononData.tt_L_out = diag(t_L*t_L');
  PhononData.tt_L_in = diag(t_L'*t_L);

  PhononData.r_L = r_L;
  PhononData.AntiXi_mode = sum(diag(r_L*r_L'));
  PhononData.rr_L_out    = diag(r_L*r_L');
  %}

  % ===== RIGHT-TO-LEFT ===
  % ===== (transmission matrix) ===
  %{
  U_minus_L = LeftLead.U_minus;                 % eigenmode matrix for out-going channels in left lead
  V_minus_L = diag(LeftLead.VecV_minus);        % eigenvelocity matrix for out-going channels in left lead
  U_plus_R_adv = RightLead.U_plus_adv;          % eigenmode matrix for in-coming channels in right lead
  V_plus_R_adv = diag(RightLead.VecV_plus_adv); % eigenvelocity matrix for in-coming channels in right lead
  t_LR = 2*1i*w*sqrt(V_minus_L)*pinv(U_minus_L)*EffG_LR* ... 
      pinv(U_plus_R_adv')*sqrt(V_plus_R_adv)/sqrt(a_R*a_L);      
  %}
  t_LR = 2*1i*w*SqrtV_minus_L*InvU_minus_L*EffG_LR* ... 
      InvTranU_plus_R_adv*SqrtV_plus_R_adv/sqrt(a_R*a_L);        
      
  % ===== (reflection matrix) ===
  %{
  U_plus_R = RightLead.U_plus;          % eigenmode matrix for reflected, out-going channels in right lead
  V_plus_R = diag(RightLead.VecV_plus); % eigenvelocity matrix for reflected out-going channels in right lead   
  EdgeG  = EffG_RR;
  EdgeG0 = RightLead.MatBulkG;
  r_RR = 2*1i*w*sqrt(V_plus_R)*pinv(U_plus_R)*(EdgeG-EdgeG0)* ...
      pinv(U_plus_R_adv')*sqrt(V_plus_R_adv)/sqrt(a_R*a_R);
  %}
  r_RR = 2*1i*w*SqrtV_plus_R*InvU_plus_R*(EffG_RR-EffG0_RR)* ...
      InvTranU_plus_R_adv*SqrtV_plus_R_adv/sqrt(a_R*a_R);
  
  % clearvars U_minus_L V_minus_L U_plus_R_adv V_plus_R_adv U_plus_R V_plus_R EdgeG EdgeG0

  % PhononData.t_LR = t_LR;
  PhononData.t_LR = sparse( t_LR );
  PhononData.xi_abs_R = diag(t_LR*t_LR');
  PhononData.xi_inc_R = diag(t_LR'*t_LR);

  % PhononData.r_RR = r_RR;
  PhononData.r_RR = sparse( r_RR );
  PhononData.xi_ref_R    = diag(r_RR*r_RR');

  clearvars t_LR r_RR;

  %{
  PhononData.t_R = t_R;
  PhononData.tt_R_out = diag(t_R*t_R');
  PhononData.tt_R_in = diag(t_R'*t_R);

  PhononData.r_R = r_R;
  PhononData.rr_R_out    = diag(r_R*r_R');
  %}
    
  % === CAROLI METHOD TO CALCULATE TOTAL TRANSMISSION ===        
  Sigma_R = HCR*SurfGR*HRC; % self-energy matrix corresponding to interaction with right lead 
  Sigma_L = HCL*SurfGL*HLC; % self-energy matrix corresponding to interaction with left lead
  Gamma_R = 1i*(Sigma_R-Sigma_R'); % right-lead dissipation matrix
  Gamma_L = 1i*(Sigma_L-Sigma_L'); % left-lead dissipation matrix
    
  TempMatrix = Gamma_L*EffG_LR*Gamma_R*EffG_LR';
  PhononData.Xi_negf = trace(TempMatrix);          

  %{
  % save('LDOS','LDOSData','wvec'); % DEBUG
  %}
end 
% end 


