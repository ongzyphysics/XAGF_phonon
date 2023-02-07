function WritePhononTransmissionOutput(DataFilesDir,wvec,LeftPhonon,RightPhonon,PhononData)
% NOTE: We write phonon and reflection coefficient data for individual modes
% NOTE: Format of "Output_Left_Transmission.dat":
% NOTE:   Frequency, band index, k_x, k_y, k_z, v_x, max t_x, t_x (transmission prob.), v_x, v_y, v_z
% NOTE: Format of "Output_Left_Reflection.dat":
% NOTE:   Frequency, band index, k_x, k_y, k_z, v_x, max r_x, r_x (transmission prob.), v_x, v_y, v_z
% NOTE: Format of "Output_Right_Transmission.dat":
% NOTE:   Frequency, band index, k_x, k_y, k_z, v_x, max t_x, t_x (transmission prob.), v_x, v_y, v_z
% NOTE: Format of "Output_Right_Reflection.dat":
% NOTE:   Frequency, band index, k_x, k_y, k_z, v_x, max r_x, r_x (transmission prob.), v_x, v_y, v_z
% NOTE: Format of "Output_Transmission.dat":
% NOTE:   Frequency, max left-lead transmission, max right-lead transmission, transmission (NEGF), transmission (sum over modes)

  CurrDir = pwd; % current directory containing AGF code .m files

  nwmax = numel(wvec); % number of frequency points

  % ------------------------------------------------------------------------------------
  % Filenames of phonon transmission and reflection data
  % ------------------------------------------------------------------------------------
  % === Filenames ===
  filename_t_L = 'Output_Left_Transmission.dat'; 
  filename_r_L = 'Output_Left_Reflection.dat';
  filename_t_R = 'Output_Right_Transmission.dat';
  filename_r_R = 'Output_Right_Reflection.dat';
  filename_tran = 'Output_Transmission.dat';

  % === File pointers to output data files ===
  cd(DataFilesDir); % change to AGF output directory
  fid_t_L = fopen(filename_t_L,'w');
  fid_r_L = fopen(filename_r_L,'w');
  fid_t_R = fopen(filename_t_R,'w');
  fid_r_R = fopen(filename_r_R,'w');
  fid_tran = fopen(filename_tran,'w');
  cd(CurrDir); % return to directory containing AGF code .m files

  % ------------------------------------------------------------------------------------
  % Write transmission and reflection data 
  % ------------------------------------------------------------------------------------
  for nw = 1:1:nwmax
      w = wvec(nw);
  
      % === Transmission data === 
      % ===== (incident/incoming left-lead phonons) ===== 
      qx_L_temp = [LeftPhonon(nw).VecQx_minus_adv]; % wave vector (x-direction) of left lead phonons  
      qy_L_temp = [LeftPhonon(nw).VecQy_minus_adv]; % wave vector (y-direction) of left lead phonons  
      qz_L_temp = [LeftPhonon(nw).VecQz_minus_adv]; % wave vector (z-direction) of left lead phonons  
      b_L_temp =  [LeftPhonon(nw).VecB_minus_adv];  % band index of left lead phonons  
      velx_L_temp = [LeftPhonon(nw).VecVelx_minus_adv]; % group velocity (x-direction) of left lead phonons  
      vely_L_temp = [LeftPhonon(nw).VecVely_minus_adv]; % group velocity (y-direction) of left lead phonons  
      velz_L_temp = [LeftPhonon(nw).VecVelz_minus_adv]; % group velocity (z-direction) of left lead phonons  

      % ===== (incident/incoming right-lead phonons) =====
      qx_R_temp = [RightPhonon(nw).VecQx_plus_adv];     % wave vector (x-direction) of right lead phonons
      qy_R_temp = [RightPhonon(nw).VecQy_plus_adv];     % wave vector (y-direction) of right lead phonons
      qz_R_temp = [RightPhonon(nw).VecQz_plus_adv];     % wave vector (z-direction) of right lead phonons
      b_R_temp =  [RightPhonon(nw).VecB_plus_adv];      % band index of right lead phonons  
      velx_R_temp = [RightPhonon(nw).VecVelx_plus_adv]; % group velocity (x-direction) of right lead phonons
      vely_R_temp = [RightPhonon(nw).VecVely_plus_adv]; % group velocity (y-direction) of right lead phonons
      velz_R_temp = [RightPhonon(nw).VecVelz_plus_adv]; % group velocity (z-direction) of right lead phonons

      q_L_temp = [LeftPhonon(nw).VecQ_minus_adv]; % wave vector (long.) of left lead phonons (b/w 0 and 2 pi)   
      q_R_temp = [RightPhonon(nw).VecQ_plus_adv]; % wave vector (long.) of right lead phonons (b/w 0 and 2 pi) 
      t_L_temp = [PhononData(nw).xi_inc_L];          % transmission coeffs. of left lead phonons (b/w 0 and 1)
      t_R_temp = [PhononData(nw).xi_inc_R];          % transmission coeffs. of right lead phonons (b/w 0 and 1)
      q_tran1_L_temp = [LeftPhonon(nw).VecQ_tran1];  % wave vector (transv. 1) of left lead phonons (b/w 0 and 2 pi) 
      q_tran2_L_temp = [LeftPhonon(nw).VecQ_tran2];  % wave vector (transv. 2) of left lead phonons (b/w 0 and 2 pi) 
      q_tran1_R_temp = [RightPhonon(nw).VecQ_tran1]; % wave vector (transv. 1) of right lead phonons (b/w 0 and 2 pi) 
      q_tran2_R_temp = [RightPhonon(nw).VecQ_tran2]; % wave vector (transv. 2) of right lead phonons (b/w 0 and 2 pi) 
      v_L_temp = [LeftPhonon(nw).VecV_minus_adv]; % group velocity of left lead phonons
      v_R_temp = [RightPhonon(nw).VecV_plus_adv]; % group velocity of right lead phonons
      p_L_temp = [LeftPhonon(nw).VecT_minus_adv]; % maximum transmission of propagating left lead phonons
      p_R_temp = [RightPhonon(nw).VecT_plus_adv]; % maximum transmission of propagating right lead phonons

      nrange_L_temp = gt(abs(b_L_temp),0); % indices for propagating incident left-lead states
      nrange_R_temp = gt(abs(b_R_temp),0); % indices for propagating incident right-lead states

      % ===== (remove non-propagating states) =====
      qx_L_temp = qx_L_temp(nrange_L_temp);
      qy_L_temp = qy_L_temp(nrange_L_temp);
      qz_L_temp = qz_L_temp(nrange_L_temp);
      b_L_temp  = b_L_temp(nrange_L_temp);
      velx_L_temp = velx_L_temp(nrange_L_temp);
      vely_L_temp = vely_L_temp(nrange_L_temp);
      velz_L_temp = velz_L_temp(nrange_L_temp);

      qx_R_temp = qx_R_temp(nrange_R_temp);
      qy_R_temp = qy_R_temp(nrange_R_temp);
      qz_R_temp = qz_R_temp(nrange_R_temp);
      b_R_temp  = b_R_temp(nrange_R_temp);
      velx_R_temp = velx_R_temp(nrange_R_temp);
      vely_R_temp = vely_R_temp(nrange_R_temp);
      velz_R_temp = velz_R_temp(nrange_R_temp);

      t_L_temp = t_L_temp(nrange_L_temp);
      t_R_temp = t_R_temp(nrange_R_temp);
      q_tran1_L_temp = q_tran1_L_temp(nrange_L_temp);
      q_tran2_L_temp = q_tran2_L_temp(nrange_L_temp);
      q_tran1_R_temp = q_tran1_R_temp(nrange_R_temp);
      q_tran2_R_temp = q_tran2_R_temp(nrange_R_temp);   
      v_L_temp = v_L_temp(nrange_L_temp);
      v_R_temp = v_R_temp(nrange_R_temp);   
      p_L_temp = p_L_temp(nrange_L_temp);
      p_R_temp = p_R_temp(nrange_R_temp);   

      q_L_temp = q_L_temp(nrange_L_temp);
      q_R_temp = q_R_temp(nrange_R_temp);

      % ===== (rationalize wave vectors by 2 pi factor for consistency) =====
      qx_L_temp = 2*pi*qx_L_temp;
      qy_L_temp = 2*pi*qy_L_temp;
      qz_L_temp = 2*pi*qz_L_temp;
      qx_R_temp = 2*pi*qx_R_temp;
      qy_R_temp = 2*pi*qy_R_temp;
      qz_R_temp = 2*pi*qz_R_temp;

      % ===== (write left-lead transmission data) =====
      if gt(numel(b_L_temp),0)
          for nq = 1:numel(b_L_temp)
              % fprintf(fid_t_L,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
              %         w, q_L_temp(nq), q_tran1_L_temp(nq), q_tran2_L_temp(nq), ... 
              %         v_L_temp(nq), p_L_temp(nq), t_L_temp(nq));            
              fprintf(fid_t_L,'%12.4e %3d %12.4e %12.4e %12.4e %12.4e %14.6e %14.6e %12.4e %12.4e %12.4e \n', ...
                      w, b_L_temp(nq), qx_L_temp(nq), qy_L_temp(nq), qz_L_temp(nq), v_L_temp(nq), ... 
                      p_L_temp(nq), t_L_temp(nq), velx_L_temp(nq), vely_L_temp(nq), velz_L_temp(nq));            
          end
      end

      % ===== (write right-lead transmission data) =====
      if gt(numel(b_R_temp),0)
          for nq = 1:numel(b_R_temp)
              % fprintf(fid_t_R,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
              %         w, q_R_temp(nq), q_tran1_R_temp(nq), q_tran2_R_temp(nq), ... 
              %         v_R_temp(nq), p_R_temp(nq), t_R_temp(nq));            
              fprintf(fid_t_R,'%12.4e %3d %12.4e %12.4e %12.4e %12.4e %14.6e %14.6e %12.4e %12.4e %12.4e \n', ...
                      w, b_R_temp(nq), qx_R_temp(nq), qy_R_temp(nq), qz_R_temp(nq), v_R_temp(nq), ... 
                      p_R_temp(nq), t_R_temp(nq), velx_R_temp(nq), vely_R_temp(nq), velz_R_temp(nq));       
          end
      end
    
      % === Reflection data ===
      % ===== (outgoing left-lead phonons) =====
      qxr_L_temp = [LeftPhonon(nw).VecQx_minus]; % wave vector (x-direction) of reflected left lead phonons  
      qyr_L_temp = [LeftPhonon(nw).VecQy_minus]; % wave vector (y-direction) of reflected left lead phonons  
      qzr_L_temp = [LeftPhonon(nw).VecQz_minus]; % wave vector (z-direction) of reflected left lead phonons  
      br_L_temp  = [LeftPhonon(nw).VecB_minus];  % band index of left lead phonons  
      velxr_L_temp = [LeftPhonon(nw).VecVelx_minus]; % group (x-direction) of reflected left lead phonons  
      velyr_L_temp = [LeftPhonon(nw).VecVely_minus]; % group (y-direction) of reflected left lead phonons  
      velzr_L_temp = [LeftPhonon(nw).VecVelz_minus]; % group (z-direction) of reflected left lead phonons  

      qr_L_temp = [LeftPhonon(nw).VecQ_minus]; % wave vector (long.) of reflected left lead phonons (b/w 0 and 2 pi)  
      tr_L_temp = [PhononData(nw).xi_ref_L];        % reflection coeffs. of left lead phonons (b/w 0 and 1)
      qr_tran1_L_temp = [LeftPhonon(nw).VecQ_tran1];  % wave vector (transv. 1) of left lead phonons (b/w 0 and 2 pi)
      qr_tran2_L_temp = [LeftPhonon(nw).VecQ_tran2];  % wave vector (transv. 2) of left lead phonons (b/w 0 and 2 pi)
      vr_L_temp = [LeftPhonon(nw).VecV_minus]; % group velocity of reflected left lead phonons
      pr_L_temp = [LeftPhonon(nw).VecT_minus]; % maximum transmission of propagating reflected left lead phonons

      % ===== (outgoing right-lead phonons) =====
      qxr_R_temp = [RightPhonon(nw).VecQx_plus]; % wave vector (x-direction) of reflected right lead phonons  
      qyr_R_temp = [RightPhonon(nw).VecQy_plus]; % wave vector (y-direction) of reflected right lead phonons  
      qzr_R_temp = [RightPhonon(nw).VecQz_plus]; % wave vector (z-direction) of reflected right lead phonons  
      br_R_temp  = [RightPhonon(nw).VecB_plus];  % band index of left lead phonons  
      velxr_R_temp = [RightPhonon(nw).VecVelx_plus]; % group (x-direction) of reflected right lead phonons  
      velyr_R_temp = [RightPhonon(nw).VecVely_plus]; % group (y-direction) of reflected right lead phonons  
      velzr_R_temp = [RightPhonon(nw).VecVelz_plus]; % group (z-direction) of reflected right lead phonons  

      qr_R_temp = [RightPhonon(nw).VecQ_plus]; % wave vector (long.) of reflected right lead phonons (b/w 0 and 2 pi)   
      tr_R_temp = [PhononData(nw).xi_ref_R];        % reflection coeffs. of right lead phonons (b/w 0 and 1)
      qr_tran1_R_temp = [RightPhonon(nw).VecQ_tran1];  % wave vector (transv. 1) of right lead phonons (b/w 0 and 2 pi)
      qr_tran2_R_temp = [RightPhonon(nw).VecQ_tran2];  % wave vector (transv. 2) of right lead phonons (b/w 0 and 2 pi)
      vr_R_temp = [RightPhonon(nw).VecV_plus]; % group velocity of reflected right lead phonons
      pr_R_temp = [RightPhonon(nw).VecT_plus]; % maximum transmission of propagating reflected right lead phonons

      nrrange_L_temp = gt(abs(br_L_temp),0); % indices for outgoing left-lead propagating states
      nrrange_R_temp = gt(abs(br_R_temp),0); % indices for outgoing right-lead propagating states

      % ===== (remove non-propagating left-lead states) =====
      qxr_L_temp = qxr_L_temp(nrrange_L_temp);
      qyr_L_temp = qyr_L_temp(nrrange_L_temp);
      qzr_L_temp = qzr_L_temp(nrrange_L_temp);
      br_L_temp  = br_L_temp(nrrange_L_temp);
      velxr_L_temp = velxr_L_temp(nrrange_L_temp);
      velyr_L_temp = velyr_L_temp(nrrange_L_temp);
      velzr_L_temp = velzr_L_temp(nrrange_L_temp);

      tr_L_temp = tr_L_temp(nrrange_L_temp);
      qr_tran1_L_temp = qr_tran1_L_temp(nrrange_L_temp);
      qr_tran2_L_temp = qr_tran2_L_temp(nrrange_L_temp);
      vr_L_temp = vr_L_temp(nrrange_L_temp);
      pr_L_temp = pr_L_temp(nrrange_L_temp);

      qr_L_temp = qr_L_temp(nrrange_L_temp);

      % ===== (remove non-propagating right-lead states) =====
      qxr_R_temp = qxr_R_temp(nrrange_R_temp);
      qyr_R_temp = qyr_R_temp(nrrange_R_temp);
      qzr_R_temp = qzr_R_temp(nrrange_R_temp);
      br_R_temp  = br_R_temp(nrrange_R_temp);
      velxr_R_temp = velxr_R_temp(nrrange_R_temp);
      velyr_R_temp = velyr_R_temp(nrrange_R_temp);
      velzr_R_temp = velzr_R_temp(nrrange_R_temp);

      tr_R_temp = tr_R_temp(nrrange_R_temp);

      qr_tran1_R_temp = qr_tran1_R_temp(nrrange_R_temp);
      qr_tran2_R_temp = qr_tran2_R_temp(nrrange_R_temp);
      vr_R_temp = vr_R_temp(nrrange_R_temp);
      pr_R_temp = pr_R_temp(nrrange_R_temp);

      qr_R_temp = qr_R_temp(nrrange_R_temp);

      % ===== (rationalize wave vectors by 2 pi factor for consistency) =====
      qxr_L_temp = 2*pi*qxr_L_temp;
      qyr_L_temp = 2*pi*qyr_L_temp;
      qzr_L_temp = 2*pi*qzr_L_temp;
      qxr_R_temp = 2*pi*qxr_R_temp;
      qyr_R_temp = 2*pi*qyr_R_temp;
      qzr_R_temp = 2*pi*qzr_R_temp;

      % ===== (write left-lead reflection data) =====
      if gt(numel(br_L_temp),0)
          for nq = 1:numel(br_L_temp)
              % fprintf(fid_r_L,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
              %         w, qr_L_temp(nq), qr_tran1_L_temp(nq), qr_tran2_L_temp(nq), ... 
              %         vr_L_temp(nq), pr_L_temp(nq), tr_L_temp(nq));            
              fprintf(fid_r_L,'%12.4e %3d %12.4e %12.4e %12.4e %12.4e %14.6e %14.6e %12.4e %12.4e %12.4e \n', ...
                      w, br_L_temp(nq), qxr_L_temp(nq), qyr_L_temp(nq), qzr_L_temp(nq), vr_L_temp(nq), ... 
                      pr_L_temp(nq), tr_L_temp(nq), velxr_L_temp(nq), velyr_L_temp(nq), velzr_L_temp(nq));            
          end
      end

      % ===== (write right-lead reflection data) =====
      if gt(numel(br_R_temp),0)
          for nq = 1:numel(br_R_temp)
              % fprintf(fid_r_R,'%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n ', ...
              %         w, qr_R_temp(nq), qr_tran1_R_temp(nq), qr_tran2_R_temp(nq), ... 
              %         vr_R_temp(nq), pr_R_temp(nq), tr_R_temp(nq));            
              fprintf(fid_r_R,'%12.4e %3d %12.4e %12.4e %12.4e %12.4e %14.6e %14.6e %12.4e %12.4e %12.4e \n', ...
                      w, br_R_temp(nq), qxr_R_temp(nq), qyr_R_temp(nq), qzr_R_temp(nq), vr_R_temp(nq), ... 
                      pr_R_temp(nq), tr_R_temp(nq), velxr_R_temp(nq), velyr_R_temp(nq), velzr_R_temp(nq));            
          end
      end

      Xi_L = real(LeftPhonon(nw).Xi_negf);
      Xi_R = real(RightPhonon(nw).Xi_negf);
      Xi_C_negf = real(PhononData(nw).Xi_negf);
      Xi_C_mode = real(PhononData(nw).Xi_mode);
    
      fprintf(fid_tran,'%16.8e %16.8e %16.8e %16.8e %16.8e \n', w, Xi_L, Xi_R, Xi_C_negf, Xi_C_mode);                
  end

  fprintf(1,'\t  <%s> \n', filename_t_L);
  fprintf(1,'\t  <%s> \n', filename_r_L);
  fprintf(1,'\t  <%s> \n', filename_t_R);
  fprintf(1,'\t  <%s> \n', filename_r_R);
  fprintf(1,'\t  <%s> \n', filename_tran);

  fclose(fid_t_L);
  fclose(fid_r_L);
  fclose(fid_t_R);
  fclose(fid_r_R);
  fclose(fid_tran);
end







