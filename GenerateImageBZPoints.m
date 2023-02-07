function [qdata_out, Udata_out, vdata_out, bdata_out] = GenerateImageBZPoints(w_in,vlong_in,qfrac1_in,qfrac2_in,qfrac3_in,p_in,Udata_in,Lead_in,LeadPhonDisp_in)
% NOTE: We implement the Boykin-Klimeck method here
% NOTE: This function takes in "Udata_in" which contains the left/right-lead eigenvectors for the AGF layer scheme
% NOTE: This function is called by "MapBulkPhononModes.m"
% NOTE: It maps the above eigenvectors to eigenvectors corresponding to bulk phonon modes of the lead
% NOTE: It also maps the non-primitive AGF wave vectors to the primitive wave vectors of the bulk phonons

  Gvec1 = Lead_in.gvec_long;  % reciprocal lattice vector (AGF, non-primitive), longitudinal direction
  Gvec2 = Lead_in.gvec_tran1; % reciprocal lattice vector (AGF, non-primitive), transverse direction 1
  Gvec3 = Lead_in.gvec_tran2; % reciprocal lattice vector (AGF, non-primitive), transverse direction 2
  PrimGvec1 = LeadPhonDisp_in.Gvec1; % reciprocal lattice vector (Brillouin Zone, primitive)
  PrimGvec2 = LeadPhonDisp_in.Gvec2; % reciprocal lattice vector (Brillouin Zone, primitive)
  PrimGvec3 = LeadPhonDisp_in.Gvec3; % reciprocal lattice vector (Brillouin Zone, primitive)

  ScaleFactor = abs(PrimGvec1*cross(PrimGvec2,PrimGvec3)')/abs(Gvec1*cross(Gvec2,Gvec3)'); 
    % ratio of AGF unit cell to primitive unit cell = zone-folding degeneracy, should be integer-value
  ScaleFactor = round(ScaleFactor); % round to integer

  nqmax = numel(p_in);                      % number of modes/channels (evanescent and propagating)
  dummy_index = 1:nqmax;                    % index of all modes 1 to nqmax
  nqrange_propagate = dummy_index(p_in==1); % index of all propagating modes (excluding evanescent modes), use for loops over prop. modes

  if eq(ScaleFactor,1) % i.e. no need for zone unfolding
      for nq = 1:1:nqmax % loop over all modes (evanescent and propagating) at frequency w
          qvec = qfrac1_in(nq)*Gvec1 + qfrac2_in(nq)*Gvec2 + qfrac3_in(nq)*Gvec3; % 'folded' wave vector of mode
          qdata_out(nq,1:3) = qvec;
      end
     
      Udata_out = Udata_in;      
      PrimCell = LeadPhonDisp_in.PrimCell;
      [vdata_out, bdata_out] = GetModeVelocity(w_in,PrimCell,Udata_out,qdata_out,nqrange_propagate);
      return; % exit function since no need to do zone unfolding
  end

  % -------------------------------------------------------------------------------
  % Find folded and unfolded wave vectors of each extended mode 
  % Folded wave vector is stored in BZPointSet(1).qdata
  % Unfolded wave vectors are stored in BZPointSet(2).qdata, BZPointSet(3).qdata 
  %   and so on, up to BZPointSet(n).qdata where n = ScaleFactor
  % Corresponding supercell reciprocal lattice vectors 'G' are stored in 
  %   BZPointSet(1).Gdata, BZPointSet(2).Gdata and so on
  % -------------------------------------------------------------------------------

  [mode_degen, BZPointSet] = SetupBZPoint(qfrac1_in,qfrac2_in,qfrac3_in,Gvec1,Gvec2,Gvec3,PrimGvec1,PrimGvec2,PrimGvec3,ScaleFactor,nqrange_propagate);
    % mode_degen stores degeneracy of each mode (evanescent and propagating) in AGF BZ
    % BZPointSet(ndegen).qdata(nq,:) stores the ndegen-th possible unfolded wave vector at qdata(nq,:)
    % BZPointSet(ndegen).Gdata(nq,:) stores the AGF reciprocal vector from folded to unfolded wave vector


  % -------------------------------------------------------------------------------
  % Find degeneracy of each phonon mode in folded AGF BZ 
  % By "degenerate", we mean modes with the same 'folded' wave vector
  % The variable "mode_main_index" stores the indices of the propagating modes
  % However, if propagating modes at nq1, nq2 and so on are degenerate and nq1<nq2<..., 
  %   then mode_main_index(nq1) = nq1
  %        mode_main_index(nq2) = nq1
  %        mode_main_index(...) = nq1
  % The list of the indices for these degenerate modes can be obtained from: 
  %   "dummy_index(mode_main_index==nq1)"
  % -------------------------------------------------------------------------------

  mode_main_index = SetupModeMainIndex(vlong_in,qfrac1_in,qfrac2_in,qfrac3_in,Gvec1,Gvec2,Gvec3,ScaleFactor,nqrange_propagate);

  PrimCell = LeadPhonDisp_in.PrimCell;
  %{
  BZPointSet = FindCorrectBZPoint(w_in,PrimCell,nqrange_propagate,BZPointSet,mode_main_index);

  n_temp = 0;
  for nq = nqrange_propagate
      if nq==mode_main_index(nq)
          for ndegen = 1:ScaleFactor
              n_temp = n_temp + BZPointSet(ndegen).nmodes(nq);
          end
      end
  end

  if numel(nqrange_propagate)~=n_temp
      error(sprintf('<!!!> No. of unfolded modes (%d) is not equal to no. of folded modes (%d).',n_temp,numel(nqrange_propagate)));
  end
  clear n_temp;
  %}


  % -------------------------------------------------------------------------------
  % Store pre-processed data before processing in next section
  % -------------------------------------------------------------------------------

  Udata_out = Udata_in;

  for nq = 1:nqmax
      qvec = qfrac1_in(nq)*Gvec1 + qfrac2_in(nq)*Gvec2 + qfrac3_in(nq)*Gvec3; % 'folded' wave vector of mode
      qdata_out(nq,1:3) = qvec;
  end



  % -------------------------------------------------------------------------------
  % Direct lattice vector of primitive unit cell within AGF unit cell where
  %   ux, uz, uz are the components of the lattice vector
  % Number of rows in boykin_rho_vector should be number of DOFs in AGF unit cell
  % -------------------------------------------------------------------------------

  ux = reshape(repmat(LeadPhonDisp_in.BasisAtomRxvec,[3 1]),[numel(repmat(LeadPhonDisp_in.BasisAtomRxvec,[3 1])) 1]);
  uy = reshape(repmat(LeadPhonDisp_in.BasisAtomRyvec,[3 1]),[numel(repmat(LeadPhonDisp_in.BasisAtomRyvec,[3 1])) 1]);
  uz = reshape(repmat(LeadPhonDisp_in.BasisAtomRzvec,[3 1]),[numel(repmat(LeadPhonDisp_in.BasisAtomRzvec,[3 1])) 1]);
  boykin_rho_vector = [ux uy uz];
  % disp(boykin_rho_vector); % DEBUG

  nsize_super = numel(ux);             % no. of DOFs in AGF unit cell
  nsize_primi = numel(ux)/ScaleFactor; % no. of DOFs in primive unit cell

  for ndegen = 1:ScaleFactor % loop over all sets of image points
      BZPointSet(ndegen).Uq_data = zeros(nsize_primi,nqmax);
  end

  % [qdata_dummy,Udata_dummy] = GetUnfoldedEigenmodes(qdata_out,Udata_out,boykin_rho_vector,nqrange_propagate,mode_main_index,ScaleFactor,BZPointSet,PrimCell,w_in);

  % -------------------------------------------------------------------------------
  % We use the Boykin-Klimeck method to extract the unfolded eigenmodes
  % -------------------------------------------------------------------------------

  [qdata_out,Udata_out] = BoykinKlimeckMethod(BZPointSet,nqrange_propagate,qfrac1_in,qfrac2_in,qfrac3_in,Gvec1,Gvec2,Gvec3,Udata_in,boykin_rho_vector,nsize_super,nsize_primi,mode_main_index); 

  % -----------------------------------------------------------------
  % We take the eigenvectors determined from the Boykin-Klimeck 
  %   step and use them to 'fill out' the real-space eigenvectors
  % -----------------------------------------------------------------

  for ncell = 1:1:numel(Lead_in.TransCell)
      rvec = Lead_in.TransCell(ncell).rvec; 
      nrange = Lead_in.TransCell(ncell).nrange;

      for nq = nqrange_propagate
          qvec = qdata_out(nq,1:3);
          phasefactor = exp(2.0*pi*1i*qvec*rvec');
          Udata_out(nrange,nq) = Udata_out(1:nsize_super,nq)*phasefactor;
      end
  end

  for nq = nqrange_propagate 
      Udata_out(:,nq) = Udata_out(:,nq)/norm(Udata_out(:,nq)); % normalization of mode eigenvector
  end

  Udata_out = AdjustEigenmodes(qdata_out,Udata_out,nqrange_propagate);

  % -----------------------------------------------------------------
  % We find the group velocity and band index of the modes
  % -----------------------------------------------------------------
  PrimCell = LeadPhonDisp_in.PrimCell;
  [vdata_out, bdata_out] = GetModeVelocity(w_in,PrimCell,Udata_out,qdata_out,nqrange_propagate);

end


% ================================================================
% ================================================================
function [qdata_out,Udata_out] = BoykinKlimeckMethod(BZPointSet,nqrange_propagate,qfrac1_in,qfrac2_in,qfrac3_in,Gvec1,Gvec2,Gvec3,Udata_in,boykin_rho_vector,nsize_super,nsize_primi,mode_main_index); 
  ScaleFactor = nsize_super/nsize_primi;

  % PART 1: BOYKIN-KLIMECK 
  for nq = nqrange_propagate % loop over propagating modes
      qvec = qfrac1_in(nq)*Gvec1 + qfrac2_in(nq)*Gvec2 + qfrac3_in(nq)*Gvec3; % 'folded' wave vector of mode, a (1 x 3) matrix
      boykin_K_vector = repmat(qvec,[nsize_super 1]); % a (nsize_super x 3) matrix
      boykin_beta_K_vector = Udata_in(1:nsize_super,nq); % eigenvector of mode in AGF supercell, a (nsize_super x 1) matrix 
      boykin_BK_vector = exp(-2.0*pi*1i*sum(boykin_rho_vector.*boykin_K_vector,2)).*boykin_beta_K_vector;
      boykin_U_matrix = zeros(nsize_super,nsize_super); % matrix for inverse Fourier transform in Boykin-Klimeck method

      U_matrix = zeros(ScaleFactor,ScaleFactor);
      for ndegen1 = 1:ScaleFactor
          rho_vector = boykin_rho_vector(1+(ndegen1-1)*nsize_primi,1:3);
          for ndegen2 = 1:ScaleFactor
              G_vector = BZPointSet(ndegen2).Gdata(nq,:);
              phasefactor = exp(2.0*pi*1i*sum(rho_vector.*G_vector));
              U_matrix(ndegen1,ndegen2) = phasefactor*1/sqrt(ScaleFactor);
          end
      end 

      CK_vector = zeros(nsize_super,1);
      for ns = 1:nsize_primi % loop over degrees of freedom in primitive unit
          nrange = ns:nsize_primi:nsize_super;
          CK_vector(nrange,1) = inv(U_matrix)*boykin_BK_vector(nrange,1);      
      end

      boykin_CK_vector = CK_vector/norm(CK_vector);

      for ndegen = 1:ScaleFactor
          BZPointSet(ndegen).Uq_data(:,nq) = boykin_CK_vector((1:nsize_primi)+(ndegen-1)*nsize_primi,1);
      end
  end

  nqmax = size(Udata_in,1);
  dummy_index = 1:nqmax;
  Udata_out = Udata_in;

  for nq = 1:nqmax
      qvec = qfrac1_in(nq)*Gvec1 + qfrac2_in(nq)*Gvec2 + qfrac3_in(nq)*Gvec3; % 'folded' wave vector of mode, a (1 x 3) matrix
      qdata_out(nq,1:3) = qvec;
  end
  
  for ndegen = 1:ScaleFactor
      BZPointSet(ndegen).nmodes = zeros(nqmax,1);
  end

  % PART 2: SORT UNFOLDED EIGENVECTORS
  for nq = nqrange_propagate % loop over propagating modes
      subspace_mode_index = dummy_index(mode_main_index==nq); 
        % list of degenerate modes corresponding to mode nq, degenerate in the sense that they share the same folded wave vector
      subspace_U_modes = [];
      subspace_q_modes = [];

      if gt(numel(subspace_mode_index),0) % i.e. subspace exists for this mode
          if gt(numel(subspace_mode_index),1) % i.e. size of subspace greater than one
              for ndegen = 1:ScaleFactor
                  A_temp = 0.0*BZPointSet(1).Uq_data(:,nq)*BZPointSet(1).Uq_data(:,nq)';
    
                  for nq1 = subspace_mode_index % loop over subspace for img point
                      Uq_temp = BZPointSet(ndegen).Uq_data(:,nq1);
                      A_temp = A_temp + Uq_temp*Uq_temp';
                  end

                  [U_temp,S_temp,V_temp] = svd(A_temp);
                  nmodes = sum(abs(diag(S_temp))>1E-3); % no. of modes at that img point

                  if gt(nmodes,0)
                      subspace_U_modes = [subspace_U_modes GramSchmidtOrthogonalize(V_temp(:,abs(diag(S_temp))>1E-3))]; % eigenvectors of modes at img point
                      subspace_q_modes = [subspace_q_modes; repmat(BZPointSet(ndegen).qdata(nq,:),[nmodes 1])];
                  end
                  BZPointSet(ndegen).nmodes(nq) = nmodes; 
              end
          elseif eq(numel(subspace_mode_index),1) % i.e. size of subspace is one
              for ndegen = 1:ScaleFactor
                  Uq_temp = BZPointSet(ndegen).Uq_data(:,nq);
                  q_temp = BZPointSet(ndegen).qdata(nq,:);

                  if gt(norm(Uq_temp),0.75)
                      subspace_U_modes = [subspace_U_modes Uq_temp];
                      subspace_q_modes = [subspace_q_modes; q_temp];
                  end
              end
          end

          if numel(subspace_mode_index)~=size(subspace_U_modes,2)
              warning('Possible error!!!');
              errmsg = sprintf('Size of degenerate-mode subspace is not equal to number of extracted eigenmodes in subspace!!!');
              error(errmsg);
          else
              % disp([nq numel(subspace_mode_index) size(subspace_U_modes,2) size(subspace_q_modes,1)]); % DEBUG
          end

          for nq1 = 1:numel(subspace_mode_index)
              qvec_temp = subspace_q_modes(nq1,:);
              Uq_temp = subspace_U_modes(:,nq1);

              boykin_q_vector = repmat(qvec_temp,[nsize_super 1]);
              phasefactor = exp(2.0*pi*1i*sum(boykin_rho_vector.*boykin_q_vector,2));
                % phase factor due to relative position of primitive unit cell within AGF supercell
              super_Uq_temp(1:nsize_super,1) = repmat(Uq_temp(1:nsize_primi,1),[ScaleFactor 1]);
              super_Uq_temp = super_Uq_temp.*phasefactor; 

              qdata_out(subspace_mode_index(nq1),1:3) = qvec_temp; % store 'unfolded' wave vector
              Udata_out(1:nsize_super,subspace_mode_index(nq1)) = super_Uq_temp; % store mode eigenvector (not normalized)
          end
      end
  end

  % PART 3: FIX PHASE OF FIRST PRIMITIVE UNIT CELL
  for nq = nqrange_propagate
      Uq_temp = Udata_out(1:nsize_primi,nq);

      for na = 1:nsize_primi
          norm_factor = 1.0; % normalization factor
          if gt(abs(Uq_temp(na,1))/norm(Uq_temp),0.1/sqrt(nsize_primi))
              norm_factor = abs(Uq_temp(na,1))/Uq_temp(na,1);
              break;
          end
      end       

      Udata_out(1:nsize_super,nq) = norm_factor*Udata_out(1:nsize_super,nq);
        % We adjust the eigenvector such that the first large eigenvector element is a positive real
  end
end

% ================================================================
% ================================================================
function Udata_out = AdjustEigenmodes(qdata_in,Udata_in,nqrange_propagate)
% This function looks for the column vectors in Udata_in and orthogonalizes 
% the vectors that share the same wave vector (in rows of qdata_in)
  Udata_out = Udata_in;
  qdata_out = qdata_in;

  nqmax = size(qdata_out,1);
  mode_main_index = zeros(1,nqmax);
  dummy_index = 1:nqmax;

  for nq1 = nqrange_propagate % loop over all propagating modes
      qvec1 = qdata_out(nq1,:);

      if gt(mode_main_index(nq1),0) % i.e. degeneracy has been determined for this mode
          continue
      else
          mode_main_index(nq1) = nq1;
      end

      for nq2 = nqrange_propagate % loop over all propagating modes
          qvec2 = qdata_out(nq2,:);

          if and( lt(norm(qvec1-qvec2)/norm(qvec1),1E-4), gt(nq2,nq1)) % i.e. modes nq1 and nq2 are degenerate and n2 > n1
              mode_main_index(nq2) = nq1; % set this mode as being one degenerate with nq2
          end
      end
  end

  for nq1 = nqrange_propagate % loop over all propagating modes
      if gt(sum(mode_main_index==nq1),1) % i.e. there are modes degenerate with nq1-mode
          Udata_mode = Udata_in(:,mode_main_index==nq1); % get the modes degenerate with nq1-mode
          Udata_out(:,mode_main_index==nq1) = GramSchmidtOrthogonalize(Udata_mode); % store the orthonormalized modes
      end
  end
end

% ================================================================
% ================================================================
function U_out = GramSchmidtOrthogonalize(U_in)
% This function orthonormalizes the column vectors in U_in
  U_out = U_in;
  n_subspace = size(U_out,2); % no. of vectors

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

% ================================================================
% ================================================================
function [qdata_out,Udata_out] = GetUnfoldedEigenmodes(qdata_in,Udata_in,boykin_rho_vector,nqrange_propagate,mode_main_index,ScaleFactor,BZPointSet,PrimCell,w_in)
  Udata_out = Udata_in;
  qdata_out = qdata_in;

  nsize_primi = size(BZPointSet(1).Uq_data,1);
  nsize_super = nsize_primi*ScaleFactor; 
  nqmax = size(BZPointSet(1).Uq_data,2);
  dummy_index = 1:nqmax;

  for nq = nqrange_propagate % loop over propagating modes
      if nq==mode_main_index(nq) % i.e. mode is the main/identifying mode of the AGF degenerate class of modes
          Uqdata_mode = []; % for storing bulk eigenmodes
          qdata_mode = [];  % for storing bulk wave vectors          

          for ndegen = 1:ScaleFactor % loop over img points in BZ
              nmodes = BZPointSet(ndegen).nmodes(nq); % degeneracy of img point
              if gt(nmodes,0) 
                  qvec = BZPointSet(ndegen).qdata(nq,:); % wave vector of img point

                  Hq = PrimCell(1).H*exp(2*pi*1i*PrimCell(1).rvec*qvec'); % dynamical matrix
                  for np = 2:numel(PrimCell)
                      Hq = Hq + PrimCell(np).H*exp(2*pi*1i*PrimCell(np).rvec*qvec');
                  end
                  [Uq,wwq] = eig(Hq); 
                  wq = sqrt(diag(wwq)); % get eigenfreqs
                  [wq_sort, ind_sort] = sort(abs(wq-w_in)); % sort eigenfreqs

                  for nm = 1:nmodes
                      Uqdata_mode = [Uqdata_mode Uq(:,ind_sort(nm))]; 
                      qdata_mode = [qdata_mode; qvec];
                  end
              end
          end

          if size(Uqdata_mode,2)==0
              error('<!!!> No degenerate modes!');
          end

          ns = 1;
          for nq_degen = dummy_index(mode_main_index==nq)              
              Udata_out(1:nsize_super,nq_degen) = repmat(Uqdata_mode(1:nsize_primi,ns),[ScaleFactor 1]);
              qdata_out(nq_degen,1:3) = qdata_mode(ns,1:3);
              ns = ns + 1;

              qvec_temp = qdata_out(nq_degen,1:3); % wave vector in primitive BZ
              Uq_temp = Udata_out(1:nsize_super,nq_degen);     % eigenvector in primitive BZ (not normalized)
              boykin_q_vector = repmat(qvec_temp,[nsize_super 1]);
              phasefactor = exp(2.0*pi*1i*sum(boykin_rho_vector.*boykin_q_vector,2));
                % phase factor due to relative position of primitive unit cell within AGF unit cell
              Uq_temp = Uq_temp.*phasefactor; 
              Udata_out(1:nsize_super,nq_degen) = Uq_temp;
          end
      end
  end
end


% ================================================================
% ================================================================
function BZPointSet = FindCorrectBZPoint(w_in,PrimCell,nqrange_propagate,BZPointSet_in,mode_main_index);
% This function finds the unfolded wave vectors here and associated mode degeneracy
  BZPointSet = BZPointSet_in;

  nqmax = size(BZPointSet(1).qdata,1); % no. modes/channels

  for nq = 1:nqmax
      for ndegen = 1:numel(BZPointSet)
          BZPointSet(ndegen).nmodes(nq) = 0; % for storing no. of primitive modes at that img point
      end
  end

  for nq = nqrange_propagate
      if mode_main_index(nq)==nq
          for ndegen = 1:numel(BZPointSet)
              qvec = BZPointSet(ndegen).qdata(nq,:);
    
              Hq = PrimCell(1).H*exp(2*pi*1i*PrimCell(1).rvec*qvec');
              for np = 2:numel(PrimCell)
                  Hq = Hq + PrimCell(np).H*exp(2*pi*1i*PrimCell(np).rvec*qvec');
              end
    
              wwq = real(eig(Hq));
    
              BZPointSet(ndegen).nmodes(nq) = sum(lt(abs(wwq-w_in*w_in)/w_in/w_in,1E-3)); % no. of bulk phonon modes at img point
          end
      end
  end  
end


% ================================================================
% ================================================================
function [mode_degen, BZPointSet] = SetupBZPoint(qfrac1_in,qfrac2_in,qfrac3_in,Gvec1,Gvec2,Gvec3,PrimGvec1,PrimGvec2,PrimGvec3,ScaleFactor,nqrange_propagate)
  nqmax = numel(qfrac1_in);
  mode_degen = ones(nqmax,1); % store degeneracy of each mode (evanescent and propagating) in AGF BZ
    % for propagating modes, mode_degen should be equal to ScaleFactor
    % for evanescent modes, mode_degen should be equal to one (original value) since they are irrelevant

  for ndegen = 1:ScaleFactor % loop over all sets of image points
      BZPointSet(ndegen).qdata = zeros(nqmax,3); % allocate data for folded and image wave vectors
      BZPointSet(ndegen).Gdata = zeros(nqmax,3); % allocate data for supercell reciprocal vectors for unfolding
  end

  for nq = 1:1:nqmax % loop over all modes (evanescent and propagating) at frequency w
      qvec = qfrac1_in(nq)*Gvec1 + qfrac2_in(nq)*Gvec2 + qfrac3_in(nq)*Gvec3; % 'folded' wave vector of mode

      for ndegen = 1:ScaleFactor
          BZPointSet(ndegen).qdata(nq,1:3) = qvec; % store wave vector
      end
  end

  nsrange = floor(-ScaleFactor/2):1:ceil(ScaleFactor/2); 

  for nq = nqrange_propagate % loop over propagating modes 
      qvec = BZPointSet(1).qdata(nq,1:3); % 'folded' wave vector of mode

      for ns1 = nsrange
          for ns2 = nsrange
              for ns3 = nsrange
                  qvec_img = qvec + ns1*Gvec1 + ns2*Gvec2 + ns3*Gvec3; 
                    % unfolded image point (wave vector)
                  qvec_img_bz = BrillouinZoneMapping(qvec_img,ScaleFactor,PrimGvec1,PrimGvec2,PrimGvec3); 
                    % unfolded image point in primitive BZ
                  if and( lt(norm(qvec_img-qvec_img_bz)/norm(qvec_img),1E-4), gt(norm(qvec_img-qvec)/norm(qvec_img),1E-4) )
                    % if unfolded point is in primitive BZ and distinct from folded wave vector
                      mode_degen(nq) = mode_degen(nq) + 1;
                      BZPointSet(mode_degen(nq)).qdata(nq,1:3) = qvec_img;
                  end
              end
          end
      end
  end

  % disp(mode_degen(mode_degen>1)); % DEBUG

  if( sum(mode_degen(nqrange_propagate)<ScaleFactor)>0 ) % i.e. if there are propagating modes that have mode_degen<ScaleFactor 
      warning('Possible error!!!');

      for nq = nqrange_propagate
          if mode_degen(nq)~=ScaleFactor
              errmsg = sprintf('Image point not computed correctly for mode %d at (%8.4e, %8.4e, %8.4e)!!!',...
                       nq, BZPointSet(1).qdata(nq,1), BZPointSet(1).qdata(nq,2), BZPointSet(1).qdata(nq,3));
              error(errmsg);
          end
      end
  end

  for ndegen = 1:ScaleFactor % loop over all sets of image points
      BZPointSet(ndegen).Gdata = BZPointSet(ndegen).qdata - BZPointSet(1).qdata; % set corresponding G vectors for image points in BZ
  end
end

% ================================================================
% ================================================================
function mode_main_index = SetupModeMainIndex(vlong_in,qfrac1_in,qfrac2_in,qfrac3_in,Gvec1,Gvec2,Gvec3,ScaleFactor,nqrange_propagate)
  nqmax = numel(qfrac1_in);
  mode_main_index = zeros(nqmax,1); % store index of degenerate mode that first appears, initialize to 0

  for nq1 = nqrange_propagate % loop over all propagating modes
      qvec1 = qfrac1_in(nq1)*Gvec1 + qfrac2_in(nq1)*Gvec2 + qfrac3_in(nq1)*Gvec3; % 'folded' wave vector of mode
      vlong1 = vlong_in(nq1);

      if gt(mode_main_index(nq1),0) % i.e. degeneracy has been determined for this mode
          continue
      else
          mode_main_index(nq1) = nq1;
      end

      for nq2 = nqrange_propagate % loop over all propagating modes
          qvec2 = qfrac1_in(nq2)*Gvec1 + qfrac2_in(nq2)*Gvec2 + qfrac3_in(nq2)*Gvec3; % 'folded' wave vector of mode
          vlong2 = vlong_in(nq2);
          
          if gt(abs((vlong2-vlong1)/vlong1),1E-2)
              continue;
          end

          if and( lt(norm(qvec1-qvec2)/norm(qvec1),1E-4), gt(nq2,nq1)) % i.e. modes nq1 and nq2 are degenerate and n2 > n1
              mode_main_index(nq2) = nq1; % set this mode as being one degenerate with nq2
          end
      end
  end
end

% ================================================================
% ================================================================
function [vdata_out, bdata_out] = GetModeVelocity(w_in,PrimCell_in,Udata_in,qdata_in,nqrange_propagate_in)
% This function computes the group velocity and band index for phonon modes
  nqmax = size(qdata_in,1);
  vdata_out = zeros(nqmax,3); % store phonon group velocity data
  bdata_out = zeros(nqmax,1); % store phonon band index data

  nsize_primi = size(PrimCell_in(1).H,1);  

  for nq = nqrange_propagate_in % loop over propagating modes
      qvec = qdata_in(nq,1:3); % wave vector of mode

      Hq = PrimCell_in(1).H*exp(2*pi*1i*PrimCell_in(1).rvec*qvec');
      for np = 2:numel(PrimCell_in)
          Hq = Hq + PrimCell_in(np).H*exp(2*pi*1i*PrimCell_in(np).rvec*qvec');
      end

      U_out = Udata_in(1:nsize_primi,nq)/norm(Udata_in(1:nsize_primi,nq));

      w_out = real(w_in);

      [Uq, wwq] = eig(Hq);
      % disp(sqrt(diag(wwq))/1E13); % DEBUG
      % disp(numel(PrimCell_in)); % DEBUG
      [ww_sort, ind_sort] = sort(diag(real(wwq)));
      Uq = Uq(:,ind_sort);
      [UU_tmp, ind_tmp] = max(abs(Uq'*U_out));
      bdata_out(nq,1) = ind_tmp;

      Hq_x = PrimCell_in(1).H*exp(2*pi*1i*PrimCell_in(1).rvec*qvec')*1i*PrimCell_in(1).rvec(1);
      Hq_y = PrimCell_in(1).H*exp(2*pi*1i*PrimCell_in(1).rvec*qvec')*1i*PrimCell_in(1).rvec(2);
      Hq_z = PrimCell_in(1).H*exp(2*pi*1i*PrimCell_in(1).rvec*qvec')*1i*PrimCell_in(1).rvec(3);
      for np = 2:numel(PrimCell_in)
          Hq_x = Hq_x + PrimCell_in(np).H*exp(2*pi*1i*PrimCell_in(np).rvec*qvec') ...
                 *1i*PrimCell_in(np).rvec(1);
          Hq_y = Hq_y + PrimCell_in(np).H*exp(2*pi*1i*PrimCell_in(np).rvec*qvec') ...
                 *1i*PrimCell_in(np).rvec(2);
          Hq_z = Hq_z + PrimCell_in(np).H*exp(2*pi*1i*PrimCell_in(np).rvec*qvec') ...
                 *1i*PrimCell_in(np).rvec(3);
      end

      vdata_out(nq,1) = real(1.0/(2.0*w_out)*U_out'*Hq_x*U_out);
      vdata_out(nq,2) = real(1.0/(2.0*w_out)*U_out'*Hq_y*U_out);
      vdata_out(nq,3) = real(1.0/(2.0*w_out)*U_out'*Hq_z*U_out);
  end
end

% ================================================================
% ================================================================
function qvec_out = BrillouinZoneMapping(qvec_in,ScaleFactor_in,PrimGvec1_in,PrimGvec2_in,PrimGvec3_in)
% This function maps the wave vector "qvec_in" to a point in the primitive BZ "qvec_out"
  qvec_s = qvec_in;

  % === Map wave vector to have positive fractional coefficients in reciprocal space ===
  fracqvec1_s = qvec_s*cross(PrimGvec2_in,PrimGvec3_in)'/(PrimGvec1_in*cross(PrimGvec2_in,PrimGvec3_in)');
  fracqvec2_s = qvec_s*cross(PrimGvec3_in,PrimGvec1_in)'/(PrimGvec2_in*cross(PrimGvec3_in,PrimGvec1_in)');
  fracqvec3_s = qvec_s*cross(PrimGvec1_in,PrimGvec2_in)'/(PrimGvec3_in*cross(PrimGvec1_in,PrimGvec2_in)');

  fracqvec1_s = mod(fracqvec1_s,1.0);
  fracqvec2_s = mod(fracqvec2_s,1.0);
  fracqvec3_s = mod(fracqvec3_s,1.0);

  qvec_s = fracqvec1_s*PrimGvec1_in + fracqvec2_s*PrimGvec2_in + fracqvec3_s*PrimGvec3_in;

  for ns1 = -1:1:0
      for ns2 = -1:1:0
          for ns3 = -1:1:0
              qvec_tmp = (fracqvec1_s+ns1)*PrimGvec1_in + (fracqvec2_s+ns2)*PrimGvec2_in ...
                         + (fracqvec3_s+ns3)*PrimGvec3_in; % shifted wave vector 
              if lt(norm(qvec_tmp),0.9999*norm(qvec_s)) 
                  qvec_s = qvec_tmp;
              end
          end
      end
  end

  qvec_out = qvec_s;
end

