function LeadPhonon = ComputeSurfaceGreensFunction(w_in,ww_in,InputParam)
% NOTE: Compute the left and right-facing surface Green's function of given lead
% NOTE: "InputParam" can be either "Left" or "Right"
% NOTE: This function requires "TransverseSubspaceTransform.m" and "CombineTransverseSubspaceComponents.m"

  w = w_in;   % input frequency
  ww = ww_in; % input frequency squared with small imaginary part
  
  HC = InputParam.MatHC;
  HL = InputParam.MatHL;
  HR = InputParam.MatHR;
  D  = eye(size(HC));
  a_long = InputParam.a_long; % lattice constant in longitudinal direction
  n_tran1 = InputParam.n_tran1; % number of Fourier components in transverse direction 1
  n_tran2 = InputParam.n_tran2; % number of Fourier components in transverse direction 2
  n_tran = n_tran1*n_tran2; % total number of transverse Fourier components

  if gt(n_tran,1) % i.e. if there are multiple transverse Fourier components or ntran > 1
      SubInputParam = TransverseSubspaceTransform(InputParam); 
        % transform large matrices (HC, HL and HR) into Fourier-component submatrices
        % this step is somewhat inefficient as it is repeated for every frequency point 

      % === Compute the surface Green's function, eigenmodes and eigenvelocities for each Fourier component ===
      % tic;
      for nt1 = 1:1:n_tran1
          for nt2 = 1:1:n_tran2
              % disp([nt1 nt2]); % DEBUG              
              SubLeadPhonon(nt1,nt2) = ComputeSurfaceGreensFunction(w,ww,SubInputParam(nt1,nt2));
              % fprintf(1,'%4d %4d : %f \n', nt1, nt2, toc);
          end
      end

      LeadPhonon = CombineTransverseSubspaceComponents(SubLeadPhonon,InputParam); 
      % fprintf(1,'DEBUG: Finally %4d %4d : %f \n', nt1, nt2, toc);
        % we rebuild the surface Green's function from the surface Green's function of the Fourier components
  else % if n_tran == 1
      TauR = -HR;
      TauL = -HL;
    
      m_max = 100; % maximum number of recursive steps in self-consistent inner loop

      % disp('Flag 1'); % DEBUG            
      % === Compute frequency-dependent surface Green's functions from ww, HL, HC and HR ===
      W = ww*D-HC; % initial surface hamiltonian   
      InvW = inv(W);
      WsR_old = W - TauR*InvW*TauL;
      WsL_old = W - TauL*InvW*TauR;
      Wb_old = W - TauR*InvW*TauL - TauL*InvW*TauR;
      TauR_old = -TauR*InvW*TauR;
      TauL_old = -TauL*InvW*TauL;
          
      natomdim = size(D,1);
      if eq(natomdim,0)
          if gt(real(W^2-4*TauL*TauR),0)
              InvWsL = (W - sqrt(abs(W^2-4*TauL*TauR)))/(2*TauL*TauR);
              InvWsR = (W - sqrt(abs(W^2-4*TauR*TauL)))/(2*TauR*TauL);
          else
              InvWsL = (W - 1i*sqrt(abs(W^2-4*TauL*TauR)))/(2*TauL*TauR);
              InvWsR = (W - 1i*sqrt(abs(W^2-4*TauR*TauL)))/(2*TauR*TauL);          
          end
          InvWb = inv(W - TauL*InvWsL*TauR - TauR*InvWsR*TauL);
      else      
          for m = 1:m_max
              InvWb_old = inv(Wb_old);
              WsR_new = WsR_old - TauR_old*InvWb_old*TauL_old;
              WsL_new = WsL_old - TauL_old*InvWb_old*TauR_old;
              Wb_new = Wb_old - TauR_old*InvWb_old*TauL_old - TauL_old*InvWb_old*TauR_old;
              TauR_new = -TauR_old*InvWb_old*TauR_old;
              TauL_new = -TauL_old*InvWb_old*TauL_old;

              epsilon = norm(Wb_new-Wb_old)/norm(Wb_new);
            
              WsR_old = WsR_new;
              WsL_old = WsL_new;
              Wb_old = Wb_new;
              TauR_old = TauR_new;
              TauL_old = TauL_new;
              
              if or(ge(m,m_max),lt(epsilon,1E-8))
                  WsR = WsR_old;
                  WsL = WsL_old;
                  InvWsR = inv(WsR_new);
                  InvWsL = inv(WsL_new);
                  InvWb = inv(Wb_new);
    
                  break;
              end
          end
      end
          
      % == Compute total transmission using Caroli formula and method === 
      Gamma_R = 1i*((TauR*InvWsR*TauL)-(TauR*InvWsR*TauL)'); 
      Gamma_L = 1i*((TauL*InvWsL*TauR)-(TauL*InvWsL*TauR)');
      Xi_negf = real(trace(Gamma_L*InvWb'*Gamma_R*InvWb)); % total transmittance, equal to number of transmitting channels 
        
      % === Compute transition amplitudes and transmission coeffs using mode analysis method ===           
      F_plus = -InvWsR*TauL;            % Bloch matrix for right-moving, right-decaying retarded states
      F_plus_adv = -InvWsR'*TauL;       % Bloch matrix for left-moving, right-decaying retarded states
      InvF_minus = (-TauL*InvWsL')';    % Inverse Bloch matrix for left-moving, left-decaying advanced states
      InvF_minus_adv = (-TauL*InvWsL)'; % Inverse Bloch matrix for right-moving, left-decaying advanced states 
  
      [U_plus,L_plus] = eig(full(F_plus));                      
        % Eigenvectors and values (phase factors) of retarded right-propagating Bloch matrix, for final states (left-to-right)
      % [U_plus_adv,L_plus_adv] = eig(F_plus_adv);          
        % Eigenvectors and values (phase factors) of advanced right-propagating Bloch matrix, for initial states (right-to-left)
      [U_minus,InvL_minus] = eig(full(InvF_minus));             
        % Eigenvectors and values (phase factors) of retarded left-propagating Bloch matrix, for final states (right-to-left)
      % [U_minus_adv,InvL_minus_adv] = eig(InvF_minus_adv); 
        % Eigenvectors and values (phase factors) of advanced left-propagating Bloch matrix, for initial states (left-to-right)
         
      tol_L = 1E-3; % tolerance for distinguishing propagating modes from evanescent modes
      Q_plus = imag(log(diag(L_plus))).* ...
               and(gt(diag(abs(L_plus)),1-tol_L),lt(diag(abs(L_plus)),1+tol_L)); % column vector
        % Extract wave vectors of right-propagating retarded states (between 0 and 2 pi), set to zero for evanescent states
      % Q_plus_adv = imag(log(diag(L_plus_adv))).* ... 
      %              and(gt(diag(abs(L_plus_adv)),1-tol_L),lt(diag(abs(L_plus_adv)),1+tol_L)); % column vector
        % Extract wave vectors of right-propagating advanced states (between 0 and 2 pi), set to zero for evanescent states
      Q_minus = -imag(log(diag(InvL_minus))).* ... 
                and(gt(diag(abs(InvL_minus)),1-tol_L),lt(diag(abs(InvL_minus)),1+tol_L)); % column vector
        % Extract wave vectors of left-propagating retarded states (between 0 and 2 pi), set to zero for evanescent states
      % Q_minus_adv = -imag(log(diag(InvL_minus_adv))).* ... 
      %               and(gt(diag(abs(InvL_minus_adv)),1-tol_L),lt(diag(abs(InvL_minus_adv)),1+tol_L)); % column vector
        % Extract wave vectors of left-propagating advanced states (between 0 and 2 pi), set to zero for evanescent states

      U_plus      = OrthonormEigenmodes(U_plus,Q_plus);
      % U_plus_adv  = OrthonormEigenmodes(U_plus_adv,Q_plus_adv);
      U_minus     = OrthonormEigenmodes(U_minus,Q_minus);
      % U_minus_adv = OrthonormEigenmodes(U_minus_adv,Q_minus_adv);

      U_plus_adv = U_plus;   % Eigenvectors of advanced right-propagating Bloch matrix, for initial states (right-to-left)
      U_minus_adv = U_minus; % Eigenvectors of advanced left-propagating Bloch matrix, for initial states (left-to-right)
      Q_plus_adv = Q_plus;   % wave vectors of right-propagating advanced states (between 0 and 2 pi) 
      Q_minus_adv = Q_minus; % wave vectors of left-propagating advanced states (between 0 and 2 pi) 

      nrange_plus_ev = le(abs(Q_plus),0); % Indices of evanescent modes (left-to-right)
      nrange_plus_pr = gt(abs(Q_plus),0); % Indices of propagating modes (left-to-right)
      nrange_minus_ev = le(abs(Q_minus),0); % Indices of evanescent modes (right-to-left) 
      nrange_minus_pr = gt(abs(Q_minus),0); % Indices of propagating modes (right-to-left) 

      U_plus_adv(:,nrange_plus_pr) = U_minus(:,nrange_minus_pr);   % Set advanced propagating states
      U_plus_adv(:,nrange_plus_ev) = U_plus(:,nrange_plus_ev);     % Set advanced evanescent states
      U_minus_adv(:,nrange_minus_pr) = U_plus(:,nrange_plus_pr);   % Set advanced propagating states
      U_minus_adv(:,nrange_minus_ev) = U_minus(:,nrange_minus_ev); % Set advanced evanescent states

      Q_plus_adv(nrange_plus_pr,:) = Q_minus(nrange_minus_pr,:);   % Set wave vectors (advanced propagating states)
      Q_plus_adv(nrange_plus_ev,:) = Q_plus(nrange_plus_ev,:);     % Set wave vectors (advanced evanescent states)
      Q_minus_adv(nrange_minus_pr,:) = Q_plus(nrange_plus_pr,:);   % Set wave vectors (advanced propagating states)
      Q_minus_adv(nrange_minus_ev,:) = Q_minus(nrange_minus_ev,:); % Set wave vectors (advanced evanescent states)

      V_plus = real(a_long/(2*w)*U_plus'*Gamma_R*U_plus); 
        % Group velocities of right-propagating retarded states, should be positive       
      V_plus_adv = -real(a_long/(2*w)*U_plus_adv'*Gamma_R*U_plus_adv);
        % Group velocities of right-propagating advanced states, should be negative
      V_minus = -real(a_long/(2*w)*U_minus'*Gamma_L*U_minus); 
        % Group velocities of left-propagating retarded states, should be negative       
      V_minus_adv = real(a_long/(2*w)*U_minus_adv'*Gamma_L*U_minus_adv);
        % Group velocities of left-propagating advanced states, should be positive

      VecV_plus = diag(V_plus);
      VecV_plus_adv = diag(V_plus_adv);
      VecV_minus = diag(V_minus); 
      VecV_minus_adv = diag(V_minus_adv); 
      VecV_plus(nrange_plus_ev,1) = 0; % set velocity of evanescent states to zero
      VecV_plus_adv(nrange_plus_ev,1) = 0; % set velocity of evanescent states to zero
      VecV_minus(nrange_minus_ev,1) = 0; % set velocity of evanescent states to zero
      VecV_minus_adv(nrange_minus_ev,1) = 0; % set velocity of evanescent states to zero
      
      V_plus = sparse( diag(VecV_plus) ); % sparsify matrix
      V_plus_adv = sparse( diag(VecV_plus_adv) ); % sparsify matrix      
      V_minus = sparse( diag( VecV_minus) ); % sparsify matrix
      V_minus_adv = sparse( diag(VecV_minus_adv) ); % sparsify matrix     
      
      %{
      BigG = InvWb;    
      t_L = 2*1i*w*sqrt(V_plus)*pinv(U_plus)*BigG*pinv(U_minus_adv')*sqrt(V_minus_adv)/a_long;
      T_plus = diag(t_L*t_L');      % Perfect left-to-right transmission (out-going) in lead 
      T_minus_adv = diag(t_L'*t_L); % Perfect left-to-right transmission (in-coming) in lead 
      t_R = 2*1i*w*sqrt(V_minus)*pinv(U_minus)*BigG*pinv(U_plus_adv')*sqrt(V_plus_adv)/a_long;
      T_minus = diag(t_R*t_R');     % Perfect right-to-left transmission (out-going) in lead 
      T_plus_adv = diag(t_R'*t_R);  % Perfect right-to-left transmission (in-coming) in lead 
      %}

      % === Store computed variables === 
      LeadPhonon.MatSurfGR = InvWsR; % Right surface Green's function
      LeadPhonon.MatSurfGL = InvWsL; % Left surface Green's function
      LeadPhonon.MatBulkG  = InvWb;  % Bulk Green's function

      LeadPhonon.U_plus      = U_plus;      % Right-propagating eigenmode matrix (retarded)
      LeadPhonon.U_plus_adv  = U_plus_adv;  % Right-propagating eigenmode matrix (advanced)
      LeadPhonon.U_minus     = U_minus;     % Left-propagating eigenmode matrix (retarded)
      LeadPhonon.U_minus_adv = U_minus_adv; % Left-propagating eigenmode matrix (advanced)

      LeadPhonon.V_plus      = V_plus;      % Velocity matrices of right-propagating eigenmodes (retarded)
      LeadPhonon.V_plus_adv  = V_plus_adv;  % Velocity matrices of right-propagating eigenmodes (advanced)
      LeadPhonon.V_minus     = V_minus;     % Velocity matrices of left-propagating eigenmodes (retarded)
      LeadPhonon.V_minus_adv = V_minus_adv; % Velocity matrices of left-propagating eigenmodes (advanced)
        
      LeadPhonon.VecV_plus      = diag(V_plus);      % 1D array of velocity of right-propagating eigenmodes (retarded)
      LeadPhonon.VecV_plus_adv  = diag(V_plus_adv);  % 1D array of velocity of right-propagating eigenmodes (advanced)
      LeadPhonon.VecV_minus     = diag(V_minus);     % 1D array of velocity of left-propagating eigenmodes (retarded)
      LeadPhonon.VecV_minus_adv = diag(V_minus_adv); % 1D array of velocity of left-propagating eigenmodes (advanced)
        
      LeadPhonon.VecQ_plus      = Q_plus;      % between 0 ans 2 pi
      LeadPhonon.VecQ_plus_adv  = Q_plus_adv;  % between 0 ans 2 pi
      LeadPhonon.VecQ_minus     = Q_minus;     % between 0 ans 2 pi
      LeadPhonon.VecQ_minus_adv = Q_minus_adv; % between 0 ans 2 pi

      LeadPhonon.VecL_plus      = diag(L_plus);      % DEBUG: Eigenvalues for retarded right-propagating Bloch matrix
      LeadPhonon.VecInvL_minus  = diag(InvL_minus);  % DEBUG: Eigenvalues for retarded left-propagating Bloch matrix

      %{
      LeadPhonon.VecT_plus      = T_plus;
      LeadPhonon.VecT_plus_adv  = T_plus_adv;
      LeadPhonon.VecT_minus     = T_minus;
      LeadPhonon.VecT_minus_adv = T_minus_adv;

      LeadPhonon.Xi_mode = real(sum(T_plus)); % number of channels using extended AGF, should be integer
      %}

      LeadPhonon.Xi_negf = Xi_negf; % number of channels using Caroli, should be integer
      LeadPhonon.VecQ_tran1 = 0.0*diag(V_plus);
      LeadPhonon.VecQ_tran2 = 0.0*diag(V_plus);
  end   

  % === Allocate variables for mapping to primitive lattice Brillouin zone === 
  LeadPhonon.VecQx_plus = [];
  LeadPhonon.VecQy_plus = [];
  LeadPhonon.VecQz_plus = [];
  LeadPhonon.VecVelx_plus = [];
  LeadPhonon.VecVely_plus = [];
  LeadPhonon.VecVelz_plus = [];
  LeadPhonon.VecDw_plus = [];
  LeadPhonon.VecB_plus  = [];

  LeadPhonon.VecQx_plus_adv = [];
  LeadPhonon.VecQy_plus_adv = [];
  LeadPhonon.VecQz_plus_adv = [];
  LeadPhonon.VecVelx_plus_adv = [];
  LeadPhonon.VecVely_plus_adv = [];
  LeadPhonon.VecVelz_plus_adv = [];
  LeadPhonon.VecDw_plus_adv = [];
  LeadPhonon.VecB_plus_adv  = [];

  LeadPhonon.VecQx_minus = [];
  LeadPhonon.VecQy_minus = [];
  LeadPhonon.VecQz_minus = [];
  LeadPhonon.VecVelx_minus = [];
  LeadPhonon.VecVely_minus = [];
  LeadPhonon.VecVelz_minus = [];
  LeadPhonon.VecDw_minus = [];
  LeadPhonon.VecB_minus  = [];

  LeadPhonon.VecQx_minus_adv = [];
  LeadPhonon.VecQy_minus_adv = [];
  LeadPhonon.VecQz_minus_adv = [];
  LeadPhonon.VecVelx_minus_adv = [];
  LeadPhonon.VecVely_minus_adv = [];
  LeadPhonon.VecVelz_minus_adv = [];
  LeadPhonon.VecDw_minus_adv = [];
  LeadPhonon.VecB_minus_adv  = [];
end 


% ================================================================
function U_out = OrthonormEigenmodes(U_in,Q_in)
% This function orthonormalizes the degenerate extended/propagating eigenmodes from the Bloch matrix
% using the Gram-Schmidt procedure

  U_out = U_in;

  nqmax = numel(Q_in); % number of modes

  for nq = 1:nqmax
      % U_out(:,nq) = U_out(:,nq)/norm(U_out(:,nq));
      if eq(Q_in(nq),0)
          % U_out(:,nq) = U_out(:,nq)*0.0; % DEBUG 
      end
  end

  FlagQ = zeros(size(Q_in));
  FlagQ(eq(abs(Q_in),0)) = 1; % set evanescent mode check flag to 1

  for nq1 = 1:nqmax      
      if eq(FlagQ(nq1),0) % i.e. mode has not been checked and is an extended mode
          FlagQ(nq1) = 1; % set check flag to 1
          q1_long = Q_in(nq1);  % longitudinal wave vector of mode 1
          U_long = U_in(:,nq1); % eigenvector of mode 1
          nqvec_long = nq1; % index of mode 1

          for nq2 = (nq1+1):nqmax % loop over remaining modes
              if eq(FlagQ(nq2),0)
                  q2_long = Q_in(nq2); % longitudinal wave vector of mode 2
                  epsilon = abs((q1_long - q2_long)/(0.5*q1_long+0.5*q2_long)); % relative diff. of long. wave vector

                  if lt(epsilon,1E-2) % i.e. mode 2 is numerically degenerate with mode 1
                      FlagQ(nq2) = 1;
                      U_long = [U_long U_in(:,nq2)]; % store eigenvector of mode 2
                      nqvec_long = [nqvec_long nq2]; % store index of mode 2
                  end
              end              
          end

          if gt(numel(nqvec_long),1)
              U_long = GramSchmidtOrthogonalize(U_long); % U_long stores the eigenvectors degenerate to mode 1
          end

          for n = 1:numel(nqvec_long)
              nq = nqvec_long(n); % index of modes in U_out that are degenerate with mode 1
              U_out(:,nq) = U_long(:,n); % replace mode-1-degenerate eigenvectors in U_out with orthonormalized eigenvectors           
          end
      end 
  end  
end


% ================================================================
function U_out = GramSchmidtOrthogonalize(U_in)
  U_out = U_in;
  n_subspace = size(U_out,2);

  for n1 = 1:n_subspace
      v = U_out(:,n1);

      for n2 = 1:(n1-1)
          u = U_out(:,n2);
          u = u/norm(u);
          v = v - (u'*v)*u;
      end

      U_out(:,n1) = v/norm(v);
  end
end


