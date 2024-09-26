function LeadPhonon = CombineTransverseSubspaceComponents(SubLeadPhonon,InputParam)
% Combine the surfacce Green's functions, eigenmodes and eigenvelocities of the subleads

  HC = InputParam.MatHC; % square matrix (lead)
  HL = InputParam.MatHL; % square matrix (lead)
  HR = InputParam.MatHR; % square matrix (lead)

  % === PRE-ALLOCATE THE VARIABLES ===
  LeadPhonon.MatSurfGR = zeros(size(HC));
  LeadPhonon.MatSurfGL = zeros(size(HC));
  LeadPhonon.MatBulkG  = zeros(size(HC));
  LeadPhonon.U_plus      = zeros(size(HC));
  LeadPhonon.U_plus_adv  = zeros(size(HC));
  LeadPhonon.U_minus     = zeros(size(HC));
  LeadPhonon.U_minus_adv = zeros(size(HC));
  LeadPhonon.V_plus      = sparse( zeros(size(HC)) );
  LeadPhonon.V_plus_adv  = sparse( zeros(size(HC)) );
  LeadPhonon.V_minus     = sparse( zeros(size(HC)) );
  LeadPhonon.V_minus_adv = sparse( zeros(size(HC)) );
  LeadPhonon.VecV_plus      = diag(zeros(size(HC)));
  LeadPhonon.VecV_plus_adv  = diag(zeros(size(HC)));
  LeadPhonon.VecV_minus     = diag(zeros(size(HC)));
  LeadPhonon.VecV_minus_adv = diag(zeros(size(HC)));
  LeadPhonon.VecQ_plus      = diag(zeros(size(HC)));
  LeadPhonon.VecQ_plus_adv  = diag(zeros(size(HC)));
  LeadPhonon.VecQ_minus     = diag(zeros(size(HC)));
  LeadPhonon.VecQ_minus_adv = diag(zeros(size(HC)));
  LeadPhonon.VecT_plus      = diag(zeros(size(HC)));
  LeadPhonon.VecT_plus_adv  = diag(zeros(size(HC)));
  LeadPhonon.VecT_minus     = diag(zeros(size(HC)));
  LeadPhonon.VecT_minus_adv = diag(zeros(size(HC)));
  LeadPhonon.Xi_negf = 0;
  LeadPhonon.Xi_mode = 0;
  LeadPhonon.VecQ_tran1 = diag(zeros(size(HC)));
  LeadPhonon.VecQ_tran2 = diag(zeros(size(HC)));

  TransCell = InputParam.TransCell; % data structure for transverse cells of lead
  n_tran = numel(TransCell);        % number of transverse cells or total number of Fourier components
  n_tran1 = InputParam.n_tran1;     % number of Fourier components (transverse 1) of lead
  n_tran2 = InputParam.n_tran2;     % number of Fourier components (transverse 2) of lead
  gvec1 = InputParam.gvec_tran1;    % reciprocal lattice vector (transverse 1) of lead
  gvec2 = InputParam.gvec_tran2;    % reciprocal lattice vector (transverse 2) of lead

  %--------------------------------------------------------------
  % Set transverse layout of Green's functions  
  %--------------------------------------------------------------
  %{  
  for n1 = 1:1:n_tran
      for n2 = 1:1:n_tran
          rvec1 = TransCell(n1).rvec; % lattice coord. for 'origin' subcell in transverse direction
          rvec2 = TransCell(n2).rvec; % lattice coord. for 'current' subcell in transverse direction
          drvec = rvec2 - rvec1; % lattice coord. displacement of 'current' subcell relative to 'origin' subcell    
          atran1(n1,n2) = round(mod(drvec*gvec1',n_tran1)); % +ve integer coeff for displacement from 'origin' to 'current' in trans. direction 1
          atran2(n1,n2) = round(mod(drvec*gvec2',n_tran2)); % +ve integer coeff for displacement from 'origin' to 'current' in trans. direction 2          
      end
  end


  for n1 = 1:1:n_tran
      for n2 = 1:1:n_tran
    	  if eq(n1,1) % first row of Greens function matrix
              n1_img(1,n2) = 1;  % image index of 'origin' subcell	  
              n2_img(1,n2) = n2; % image index of 'current' subcell
     	  else
     	      n1_img(n1,n2) = n1_img(1,and(atran1(1,1:n_tran)==atran1(n1,n2),atran2(1,1:n_tran)==atran2(n1,n2)));
     	      n2_img(n1,n2) = n2_img(1,and(atran1(1,1:n_tran)==atran1(n1,n2),atran2(1,1:n_tran)==atran2(n1,n2)));     	      
    	  end
      end
  end
  %}

  n1_img = InputParam.n1_img;
  n2_img = InputParam.n2_img;  
    
  %--------------------------------------------------------------
  % We obtain the surface and bulk Green's function here after 
  % the inverse Fourier transform for the inverse direction
  %--------------------------------------------------------------

  for n1 = 1:1:n_tran % loop over transverse cell
      if eq(n1,1)
          for n2 = 1:1:n_tran % loop over transverse cell
              nrange1 = TransCell(n1).nrange; % matrix indices for 'origin' subcell in transverse direction
              nrange2 = TransCell(n2).nrange; % matrix indices for 'current' subcell in transverse direction

              rvec1 = TransCell(n1).rvec; % lattice coord. for 'origin' subcell in transverse direction
              rvec2 = TransCell(n2).rvec; % lattice coord. for 'current' subcell in transverse direction
              drvec = rvec2 - rvec1; % lattice coord. displacement of 'current' subcell relative to 'origin' subcell    


              for nt1 = 1:1:n_tran1 % loop over Fourier component (transverse 1)
                  for nt2 = 1:1:n_tran2 % loop over Fourier component (transverse 2)
                      kvec = double(nt1)/double(n_tran1)*gvec1 + double(nt2)/double(n_tran2)*gvec2;
                      phasefactor = -2.0*pi*kvec*drvec'; % minus sign here becos we are doing an inverse Fourier transform
                      LeadPhonon.MatSurfGR(nrange1,nrange2) = LeadPhonon.MatSurfGR(nrange1,nrange2) ... 
                          + SubLeadPhonon(nt1,nt2).MatSurfGR * exp(1i*phasefactor);
                      LeadPhonon.MatSurfGL(nrange1,nrange2) = LeadPhonon.MatSurfGL(nrange1,nrange2) ... 
                          + SubLeadPhonon(nt1,nt2).MatSurfGL * exp(1i*phasefactor);
                      LeadPhonon.MatBulkG(nrange1,nrange2) = LeadPhonon.MatBulkG(nrange1,nrange2) ... 
                          + SubLeadPhonon(nt1,nt2).MatBulkG * exp(1i*phasefactor);
                  end
              end
          end
      else
          nrange1 = TransCell(n1).nrange; % matrix indices for 'origin' subcell in transverse direction
          nrange2 = [TransCell(1:1:n_tran).nrange]; % matrix indices for 'current' subcell in transverse direction
          nrange1_img = TransCell(n1_img(n1,1)).nrange;
          nrange2_img = [TransCell(n2_img(n1,1:1:n_tran)).nrange];
      
          LeadPhonon.MatSurfGR(nrange1,nrange2) = LeadPhonon.MatSurfGR(nrange1_img,nrange2_img);              
          LeadPhonon.MatSurfGL(nrange1,nrange2) = LeadPhonon.MatSurfGL(nrange1_img,nrange2_img);
          LeadPhonon.MatBulkG(nrange1,nrange2) = LeadPhonon.MatBulkG(nrange1_img,nrange2_img);      
      end
  end
  %{  
  %}
  
  %{
  for n1 = 1:1:n_tran % loop over transverse cell
      for n2 = 1:1:n_tran % loop over transverse cell
          nrange1 = TransCell(n1).nrange; % matrix indices for 'origin' subcell in transverse direction
          nrange2 = TransCell(n2).nrange; % matrix indices for 'current' subcell in transverse direction

          rvec1 = TransCell(n1).rvec; % lattice coord. for 'origin' subcell in transverse direction
          rvec2 = TransCell(n2).rvec; % lattice coord. for 'current' subcell in transverse direction
          drvec = rvec2 - rvec1; % lattice coord. displacement of 'current' subcell relative to 'origin' subcell    

          %{
          % --- Old code for reconstruction of real-space Greens function ---
          for nt1 = 1:1:n_tran1 % loop over Fourier component (transverse 1)
              for nt2 = 1:1:n_tran2 % loop over Fourier component (transverse 2)
                  kvec = double(nt1)/double(n_tran1)*gvec1 + double(nt2)/double(n_tran2)*gvec2;
                  phasefactor = -2.0*pi*kvec*drvec'; % minus sign here becos we are doing an inverse Fourier transform
                  LeadPhonon.MatSurfGR(nrange1,nrange2) = LeadPhonon.MatSurfGR(nrange1,nrange2) ... 
                      + SubLeadPhonon(nt1,nt2).MatSurfGR * exp(1i*phasefactor);
                  LeadPhonon.MatSurfGL(nrange1,nrange2) = LeadPhonon.MatSurfGL(nrange1,nrange2) ... 
                      + SubLeadPhonon(nt1,nt2).MatSurfGL * exp(1i*phasefactor);
                  LeadPhonon.MatBulkG(nrange1,nrange2) = LeadPhonon.MatBulkG(nrange1,nrange2) ... 
                      + SubLeadPhonon(nt1,nt2).MatBulkG * exp(1i*phasefactor);
              end
          end
          %}
          
          % --- New code to optimize reconstruction of real-space Greens function ---
          if eq(n1,1)
              for nt1 = 1:1:n_tran1 % loop over Fourier component (transverse 1)
                  for nt2 = 1:1:n_tran2 % loop over Fourier component (transverse 2)
                      kvec = double(nt1)/double(n_tran1)*gvec1 + double(nt2)/double(n_tran2)*gvec2;
                      phasefactor = -2.0*pi*kvec*drvec'; % minus sign here becos we are doing an inverse Fourier transform
                      LeadPhonon.MatSurfGR(nrange1,nrange2) = LeadPhonon.MatSurfGR(nrange1,nrange2) ... 
                          + SubLeadPhonon(nt1,nt2).MatSurfGR * exp(1i*phasefactor);
                      LeadPhonon.MatSurfGL(nrange1,nrange2) = LeadPhonon.MatSurfGL(nrange1,nrange2) ... 
                          + SubLeadPhonon(nt1,nt2).MatSurfGL * exp(1i*phasefactor);
                      LeadPhonon.MatBulkG(nrange1,nrange2) = LeadPhonon.MatBulkG(nrange1,nrange2) ... 
                          + SubLeadPhonon(nt1,nt2).MatBulkG * exp(1i*phasefactor);
                  end
              end
          else
              nrange1_img = TransCell(n1_img(n1,n2)).nrange;
              nrange2_img = TransCell(n2_img(n1,n2)).nrange;
              
              LeadPhonon.MatSurfGR(nrange1,nrange2) = LeadPhonon.MatSurfGR(nrange1_img,nrange2_img);              
              LeadPhonon.MatSurfGL(nrange1,nrange2) = LeadPhonon.MatSurfGL(nrange1_img,nrange2_img);
              LeadPhonon.MatBulkG(nrange1,nrange2) = LeadPhonon.MatBulkG(nrange1_img,nrange2_img);
          end
      end
  end
  %}  
  
  LeadPhonon.MatSurfGR = 1.0/double(n_tran)*LeadPhonon.MatSurfGR;
  LeadPhonon.MatSurfGL = 1.0/double(n_tran)*LeadPhonon.MatSurfGL;
  LeadPhonon.MatBulkG  = 1.0/double(n_tran)*LeadPhonon.MatBulkG;

  %--------------------------------------------------------------
  % We obtain the 'eigenvectors' 
  %--------------------------------------------------------------
  for n = 1:1:n_tran % loop over transverse cell
      rvec = TransCell(n).rvec;
      nrange = TransCell(n).nrange;
    
      U_plus = [];
      U_plus_adv = [];
      U_minus = [];
      U_minus_adv = [];

      for nt2 = 1:1:n_tran2 % loop over Fourier component (transverse 2)
          for nt1 = 1:1:n_tran1 % loop over Fourier component (transverse 1)
              kvec = double(nt1)/double(n_tran1)*gvec1 + double(nt2)/double(n_tran2)*gvec2;
              phasefactor = 2.0*pi*kvec*rvec';
              U_plus = [U_plus SubLeadPhonon(nt1,nt2).U_plus*exp(1i*phasefactor)];
              U_plus_adv = [U_plus_adv SubLeadPhonon(nt1,nt2).U_plus_adv*exp(1i*phasefactor)];
              U_minus = [U_minus SubLeadPhonon(nt1,nt2).U_minus*exp(1i*phasefactor)];
              U_minus_adv = [U_minus_adv SubLeadPhonon(nt1,nt2).U_minus_adv*exp(1i*phasefactor)];
          end
      end

      LeadPhonon.U_plus(nrange,:)      = 1.0/sqrt(double(n_tran))*U_plus;
      LeadPhonon.U_plus_adv(nrange,:)  = 1.0/sqrt(double(n_tran))*U_plus_adv;
      LeadPhonon.U_minus(nrange,:)     = 1.0/sqrt(double(n_tran))*U_minus;
      LeadPhonon.U_minus_adv(nrange,:) = 1.0/sqrt(double(n_tran))*U_minus_adv;
  end
  
  VecV_plus      = [];
  VecV_plus_adv  = [];
  VecV_minus     = [];
  VecV_minus_adv = [];

  VecQ_plus      = [];
  VecQ_plus_adv  = [];
  VecQ_minus     = [];
  VecQ_minus_adv = [];

  VecT_plus      = [];
  VecT_plus_adv  = [];
  VecT_minus     = [];
  VecT_minus_adv = [];

  Xi_negf = 0;
  Xi_mode = 0;

  VecQ_tran1 = [];
  VecQ_tran2 = [];

  Id_temp = ones(size(SubLeadPhonon(1,1).VecQ_plus));

  for nt2 = 1:1:n_tran2 % loop over Fourier component (transverse 2)
      for nt1 = 1:1:n_tran1 % loop over Fourier component (transverse 1)
          VecV_plus      = [VecV_plus; SubLeadPhonon(nt1,nt2).VecV_plus];
          VecV_plus_adv  = [VecV_plus_adv; SubLeadPhonon(nt1,nt2).VecV_plus_adv];
          VecV_minus     = [VecV_minus; SubLeadPhonon(nt1,nt2).VecV_minus];
          VecV_minus_adv = [VecV_minus_adv; SubLeadPhonon(nt1,nt2).VecV_minus_adv];

          VecQ_plus      = [VecQ_plus; SubLeadPhonon(nt1,nt2).VecQ_plus];
          VecQ_plus_adv  = [VecQ_plus_adv; SubLeadPhonon(nt1,nt2).VecQ_plus_adv];
          VecQ_minus     = [VecQ_minus; SubLeadPhonon(nt1,nt2).VecQ_minus];
          VecQ_minus_adv = [VecQ_minus_adv; SubLeadPhonon(nt1,nt2).VecQ_minus_adv];

          % VecT_plus      = [VecT_plus; SubLeadPhonon(nt1,nt2).VecT_plus];
          % VecT_plus_adv  = [VecT_plus_adv; SubLeadPhonon(nt1,nt2).VecT_plus_adv];
          % VecT_minus     = [VecT_minus; SubLeadPhonon(nt1,nt2).VecT_minus];
          % VecT_minus_adv = [VecT_minus_adv; SubLeadPhonon(nt1,nt2).VecT_minus_adv];

          Xi_negf = Xi_negf + SubLeadPhonon(nt1,nt2).Xi_negf;
          % Xi_mode = Xi_mode + SubLeadPhonon(nt1,nt2).Xi_mode;

          VecQ_tran1 = [VecQ_tran1; 2.0*pi*double(nt1)/double(n_tran1)*Id_temp];
          VecQ_tran2 = [VecQ_tran2; 2.0*pi*double(nt2)/double(n_tran2)*Id_temp];
      end
  end
  
  LeadPhonon.V_plus      = sparse( diag(VecV_plus) );
  LeadPhonon.V_plus_adv  = sparse( diag(VecV_plus_adv) );
  LeadPhonon.V_minus     = sparse( diag(VecV_minus) );
  LeadPhonon.V_minus_adv = sparse( diag(VecV_minus_adv) );

  LeadPhonon.VecV_plus      = VecV_plus;
  LeadPhonon.VecV_plus_adv  = VecV_plus_adv;
  LeadPhonon.VecV_minus     = VecV_minus;
  LeadPhonon.VecV_minus_adv = VecV_minus_adv;

  LeadPhonon.VecQ_plus      = VecQ_plus;
  LeadPhonon.VecQ_plus_adv  = VecQ_plus_adv;
  LeadPhonon.VecQ_minus     = VecQ_minus;
  LeadPhonon.VecQ_minus_adv = VecQ_minus_adv;

  % LeadPhonon.VecT_plus      = VecT_plus;
  % LeadPhonon.VecT_plus_adv  = VecT_plus_adv;
  % LeadPhonon.VecT_minus     = VecT_minus;
  % LeadPhonon.VecT_minus_adv = VecT_minus_adv;

  LeadPhonon.Xi_negf = Xi_negf;
  % LeadPhonon.Xi_mode = Xi_mode;

  LeadPhonon.VecQ_tran1 = VecQ_tran1;
  LeadPhonon.VecQ_tran2 = VecQ_tran2;
end



