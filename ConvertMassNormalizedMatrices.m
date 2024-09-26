function [Left, Center, Right] = ConvertMassNormalizedMatrices(LeftLeadParam,CenterParam,RightLeadParam)
% NOTE: Set the mass-normalized force constant matrices for left, center and right regions
% NOTE: 

  Lyr = CenterParam.Lyr; % structure containing IFC's and atomic masses in scattering region
  nmax = length(Lyr); % no. of layers in scattering region

  for n = 1:1:nmax % we set up the inverse square root mass matrix for mass-normalizing the force constant matrices 
      if isdiag(Lyr(n).MatM) % if diagonal matrix
          Lyr(n).InvSqrtM = sparse( diag(1./sqrt(diag(Lyr(n).MatM))) );
      else
          [V,D] = eig(full([Lyr(n).MatM]));
          Lyr(n).InvSqrtM = sparse( V*diag(1./sqrt(diag(D)))*V' );        
      end
  end

  % --------------------------------------------------
  % Set left lead parameters and variables
  % --------------------------------------------------  
  if isdiag(LeftLeadParam.MatM)
      Left.InvSqrtM = diag(1./sqrt(diag(LeftLeadParam.MatM)));
  else
      [V,D] = eig(full([LeftLeadParam.MatM]));
      Left.InvSqrtM = sparse( V*diag(1./sqrt(diag(D)))*V' );
  end
  Left.MatHL = Left.InvSqrtM * LeftLeadParam.MatKL * Left.InvSqrtM;
  Left.MatHC = Left.InvSqrtM * LeftLeadParam.MatKC * Left.InvSqrtM;
  Left.MatHR = Left.InvSqrtM * LeftLeadParam.MatKR * Left.InvSqrtM;

  Left.TransCell = LeftLeadParam.TransCell;
  Left.a_long = LeftLeadParam.a_long; % longitudinal lattice spacing
  Left.n_tran1 = LeftLeadParam.n_tran1; % number of Fourier components (transverse 1)
  Left.n_tran2 = LeftLeadParam.n_tran2; % number of Fourier components (transverse 2)
  Left.rvec_long = LeftLeadParam.rvec_long;   % lattice vector (longitudinal)
  Left.rvec_tran1 = LeftLeadParam.rvec_tran1; % lattice vector (transverse 1)
  Left.rvec_tran2 = LeftLeadParam.rvec_tran2; % lattice vector (transverse 2)
  Left.gvec_long  = cross(Left.rvec_tran1,Left.rvec_tran2) ... % reciprocal lattice vector (longitudinal)
                    /(Left.rvec_long*cross(Left.rvec_tran1,Left.rvec_tran2)');
  Left.gvec_tran1 = cross(Left.rvec_tran2,Left.rvec_long) ...  % reciprocal lattice vector (transverse 1)
                    /(Left.rvec_tran1*cross(Left.rvec_tran2,Left.rvec_long)');
  Left.gvec_tran2 = cross(Left.rvec_long,Left.rvec_tran1) ...  % reciprocal lattice vector (transverse 2)
                    /(Left.rvec_tran2*cross(Left.rvec_long,Left.rvec_tran1)');

  %--------------------------------------------------------------
  % Set transverse layout of left-lead Green's functions  
  %--------------------------------------------------------------    
  n_tran = Left.n_tran1*Left.n_tran2; % total no. of transverse cells
  
  for n1 = 1:1:n_tran
      for n2 = 1:1:n_tran
          rvec1 = Left.TransCell(n1).rvec; % lattice coord. for 'origin' subcell in transverse direction
          rvec2 = Left.TransCell(n2).rvec; % lattice coord. for 'current' subcell in transverse direction
          drvec = rvec2 - rvec1; % lattice coord. displacement of 'current' subcell relative to 'origin' subcell    
          atran1(n1,n2) = round(mod(drvec*Left.gvec_tran1',Left.n_tran1)); % +ve integer coeff for displacement from 'origin' to 'current' in trans. direction 1
          atran2(n1,n2) = round(mod(drvec*Left.gvec_tran2',Left.n_tran2)); % +ve integer coeff for displacement from 'origin' to 'current' in trans. direction 2          
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

  Left.n1_img = round(n1_img);
  Left.n2_img = round(n2_img);  

  % --------------------------------------------------
  % Set right lead parameters and variables
  % --------------------------------------------------  
  if isdiag(RightLeadParam.MatM)
      Right.InvSqrtM = diag(1./sqrt(diag(RightLeadParam.MatM)));
  else
      [V,D] = eig(full([RightLeadParam.MatM]));
      Right.InvSqrtM = sparse( V*diag(1./sqrt(diag(D)))*V' );
  end
  Right.MatHL = Right.InvSqrtM * RightLeadParam.MatKL * Right.InvSqrtM;
  Right.MatHC = Right.InvSqrtM * RightLeadParam.MatKC * Right.InvSqrtM;
  Right.MatHR = Right.InvSqrtM * RightLeadParam.MatKR * Right.InvSqrtM;

  Right.TransCell = RightLeadParam.TransCell;
  Right.a_long = RightLeadParam.a_long; % longitudinal lattice spacing
  Right.n_tran1 = RightLeadParam.n_tran1; % number of Fourier components (transverse 1)
  Right.n_tran2 = RightLeadParam.n_tran2; % number of Fourier components (transverse 2)
  Right.rvec_long = RightLeadParam.rvec_long;   % lattice vector (longitudinal)
  Right.rvec_tran1 = RightLeadParam.rvec_tran1; % lattice vector (transverse 1)
  Right.rvec_tran2 = RightLeadParam.rvec_tran2; % lattice vector (transverse 2)
  Right.gvec_long  = cross(Right.rvec_tran1,Right.rvec_tran2) ... % reciprocal lattice vector (longitudinal)
                     /(Right.rvec_long*cross(Right.rvec_tran1,Right.rvec_tran2)');
  Right.gvec_tran1 = cross(Right.rvec_tran2,Right.rvec_long) ...  % reciprocal lattice vector (transverse 1) 
                     /(Right.rvec_tran1*cross(Right.rvec_tran2,Right.rvec_long)');
  Right.gvec_tran2 = cross(Right.rvec_long,Right.rvec_tran1) ...  % reciprocal lattice vector (transverse 2) 
                     /(Right.rvec_tran2*cross(Right.rvec_long,Right.rvec_tran1)');
  
  %--------------------------------------------------------------
  % Set transverse layout of right-lead Green's functions  
  %--------------------------------------------------------------    
  n_tran = Right.n_tran1*Right.n_tran2; % total no. of transverse cells
  
  for n1 = 1:1:n_tran
      for n2 = 1:1:n_tran
          rvec1 = Right.TransCell(n1).rvec; % lattice coord. for 'origin' subcell in transverse direction
          rvec2 = Right.TransCell(n2).rvec; % lattice coord. for 'current' subcell in transverse direction
          drvec = rvec2 - rvec1; % lattice coord. displacement of 'current' subcell relative to 'origin' subcell    
          atran1(n1,n2) = round(mod(drvec*Right.gvec_tran1',Right.n_tran1)); % +ve integer coeff for displacement from 'origin' to 'current' in trans. direction 1
          atran2(n1,n2) = round(mod(drvec*Right.gvec_tran2',Right.n_tran2)); % +ve integer coeff for displacement from 'origin' to 'current' in trans. direction 2          
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

  Right.n1_img = round(n1_img);
  Right.n2_img = round(n2_img);  
  
  % --------------------------------------------------
  % Set center layer parameters and variables
  % --------------------------------------------------  
  if eq(nmax,1)
      Lyr(1).MatHC = Lyr(1).InvSqrtM * Lyr(1).MatKC * Lyr(1).InvSqrtM;
  else
      for n = 1:1:nmax
          Lyr(n).MatHC = Lyr(n).InvSqrtM * Lyr(n).MatKC * Lyr(n).InvSqrtM;
          if eq(n,1)
              Lyr(n).MatHR = Lyr(n).InvSqrtM * Lyr(n).MatKR * Lyr(n+1).InvSqrtM;
          elseif eq(n,nmax)
              Lyr(n).MatHL = Lyr(n).InvSqrtM * Lyr(n).MatKL * Lyr(n-1).InvSqrtM;
          else
              Lyr(n).MatHL = Lyr(n).InvSqrtM * Lyr(n).MatKL * Lyr(n-1).InvSqrtM;
              Lyr(n).MatHR = Lyr(n).InvSqrtM * Lyr(n).MatKR * Lyr(n+1).InvSqrtM;
          end
      end
  end

  HCL = Lyr(1).InvSqrtM * CenterParam.MatKCL * Left.InvSqrtM;
  HLC = HCL';
  HCR = Lyr(nmax).InvSqrtM * CenterParam.MatKCR * Right.InvSqrtM;
  HRC = HCR';

  Lyr = rmfield(Lyr,'MatM');
  Lyr = rmfield(Lyr,'InvSqrtM');
  Lyr = rmfield(Lyr,'MatKC');
  Lyr = rmfield(Lyr,'MatKL');
  Lyr = rmfield(Lyr,'MatKR');
  Left  = rmfield(Left,'InvSqrtM');
  Right = rmfield(Right,'InvSqrtM');
  
  Left.MatHL = sparse(Left.MatHL);
  Left.MatHC = sparse(Left.MatHC);
  Left.MatHR = sparse(Left.MatHR);

  Right.MatHL = sparse(Right.MatHL);
  Right.MatHC = sparse(Right.MatHC);
  Right.MatHR = sparse(Right.MatHR);
  
  for n = 1:1:nmax
      Lyr(n).MatHL = sparse(Lyr(n).MatHL);
      Lyr(n).MatHC = sparse(Lyr(n).MatHC);
      Lyr(n).MatHR = sparse(Lyr(n).MatHR);
  end
  
  HCL = sparse(HCL);
  HCR = sparse(HCR);
  HLC = sparse(HLC);
  HRC = sparse(HRC);  

  Center.Lyr = Lyr;
  Center.HCL = HCL;
  Center.HCR = HCR;
  Center.HLC = HLC;
  Center.HRC = HRC;

end

