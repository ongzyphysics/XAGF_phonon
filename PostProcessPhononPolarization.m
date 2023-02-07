function PostProcessPhononPolarization(InputFilesDir,OutputFilesDir,LeftCutoff_in,RightCutoff_in)


  hbar = 1.055E-34;
  meV = 1.602E-22;
  epsilon_left = 2.00; % cutoff parameter for relative distance from gamma point, between 0.0 and 1.0, default is 0.10
  epsilon_right = 2.00; % cutoff parameter for relative distance from gamma point, between 0.0 and 1.0, default is 0.10
  delta_w = 1E-3; % sensitivity parameter (dimensionless) for w degeneracy
  delta_q = 1E-4; % sensitivity parameter (dimensionless) for q degeneracy
  delta_v = 1E-3; % sensitivity parameter (dimensionless) for phonon group velocity

  [LeftPhonDisp, RightPhonDisp] = ReadPhononDispersionParameters(InputFilesDir);
  % nchmax = 31;
  % epsilon = 0.0;
  % nqmax = 1;
  CurrDir = pwd;

  % === Real-space lattice vector (left and right) ===
  LeftPhon.Tvec1 = LeftPhonDisp.Rvec1;
  LeftPhon.Tvec2 = LeftPhonDisp.Rvec2;
  LeftPhon.Tvec3 = LeftPhonDisp.Rvec3;
  LeftPhon.Tvec1 = LeftPhon.Tvec1';
  LeftPhon.Tvec2 = LeftPhon.Tvec2';
  LeftPhon.Tvec3 = LeftPhon.Tvec3';

  RightPhon.Tvec1 = RightPhonDisp.Rvec1;
  RightPhon.Tvec2 = RightPhonDisp.Rvec2;
  RightPhon.Tvec3 = RightPhonDisp.Rvec3;
  RightPhon.Tvec1 = RightPhon.Tvec1';
  RightPhon.Tvec2 = RightPhon.Tvec2';
  RightPhon.Tvec3 = RightPhon.Tvec3';

  % === Reciprocal-space lattice vector (left and right) ===
  LeftPhon.Qvec1 = 2*pi*cross(LeftPhon.Tvec2,LeftPhon.Tvec3)/(cross(LeftPhon.Tvec2,LeftPhon.Tvec3)'*LeftPhon.Tvec1);
  LeftPhon.Qvec2 = 2*pi*cross(LeftPhon.Tvec3,LeftPhon.Tvec1)/(cross(LeftPhon.Tvec3,LeftPhon.Tvec1)'*LeftPhon.Tvec2);
  LeftPhon.Qvec3 = 2*pi*cross(LeftPhon.Tvec1,LeftPhon.Tvec2)/(cross(LeftPhon.Tvec1,LeftPhon.Tvec2)'*LeftPhon.Tvec3);
  RightPhon.Qvec1 = 2*pi*cross(RightPhon.Tvec2,RightPhon.Tvec3)/(cross(RightPhon.Tvec2,RightPhon.Tvec3)'*RightPhon.Tvec1);
  RightPhon.Qvec2 = 2*pi*cross(RightPhon.Tvec3,RightPhon.Tvec1)/(cross(RightPhon.Tvec3,RightPhon.Tvec1)'*RightPhon.Tvec2);
  RightPhon.Qvec3 = 2*pi*cross(RightPhon.Tvec1,RightPhon.Tvec2)/(cross(RightPhon.Tvec1,RightPhon.Tvec2)'*RightPhon.Tvec3);

  LeftPhon.Qvec1 = 2*pi*LeftPhonDisp.Gvec1; 
  LeftPhon.Qvec2 = 2*pi*LeftPhonDisp.Gvec2;
  LeftPhon.Qvec3 = 2*pi*LeftPhonDisp.Gvec3;
  RightPhon.Qvec1 = 2*pi*RightPhonDisp.Gvec1;
  RightPhon.Qvec2 = 2*pi*RightPhonDisp.Gvec2;
  RightPhon.Qvec3 = 2*pi*RightPhonDisp.Gvec3;

  % === Shell IFC matrices (right) ===
  for n = 1:numel(LeftPhonDisp.PrimCell)
      LeftPhon.S(n).H = LeftPhonDisp.PrimCell(n).H;
      LeftPhon.S(n).rvec = LeftPhonDisp.PrimCell(n).rvec';
  end


  for n = 1:numel(RightPhonDisp.PrimCell)
      RightPhon.S(n).H = RightPhonDisp.PrimCell(n).H;
      RightPhon.S(n).rvec = RightPhonDisp.PrimCell(n).rvec';
  end

  % === Determine cutoffs (left and right) ===
  % LeftPhon.qcutoff_G  = epsilon_left*mean([norm(LeftPhon.Qvec1) norm(LeftPhon.Qvec2) norm(LeftPhon.Qvec3)]);
  % RightPhon.qcutoff_G = epsilon_right*mean([norm(RightPhon.Qvec1) norm(RightPhon.Qvec2) norm(RightPhon.Qvec3)]);
  LeftPhon.qcutoff_G  = LeftCutoff_in;
  RightPhon.qcutoff_G = RightCutoff_in;

  fprintf(1,'\n');
  fprintf(1,'++ Left cutoff wave vector is %12.8e m-1 \n',LeftPhon.qcutoff_G);
  fprintf(1,'\t G_1 = (%12.8e,%12.8e,%12.8e), length = %12.8e \n',LeftPhon.Qvec1(1),LeftPhon.Qvec1(2),LeftPhon.Qvec1(3),norm(LeftPhon.Qvec1));
  fprintf(1,'\t G_2 = (%12.8e,%12.8e,%12.8e), length = %12.8e \n',LeftPhon.Qvec2(1),LeftPhon.Qvec2(2),LeftPhon.Qvec2(3),norm(LeftPhon.Qvec2));
  fprintf(1,'\t G_3 = (%12.8e,%12.8e,%12.8e), length = %12.8e \n',LeftPhon.Qvec3(1),LeftPhon.Qvec3(2),LeftPhon.Qvec3(3),norm(LeftPhon.Qvec3));
  fprintf(1,'++ Right cutoff wave vector is %12.8e m-1 \n',RightPhon.qcutoff_G);
  fprintf(1,'\t G_1 = (%12.8e,%12.8e,%12.8e), length = %12.8e \n',RightPhon.Qvec1(1),RightPhon.Qvec1(2),RightPhon.Qvec1(3),norm(RightPhon.Qvec1));
  fprintf(1,'\t G_2 = (%12.8e,%12.8e,%12.8e), length = %12.8e \n',RightPhon.Qvec2(1),RightPhon.Qvec2(2),RightPhon.Qvec2(3),norm(RightPhon.Qvec2));
  fprintf(1,'\t G_3 = (%12.8e,%12.8e,%12.8e), length = %12.8e \n',RightPhon.Qvec3(1),RightPhon.Qvec3(2),RightPhon.Qvec3(3),norm(RightPhon.Qvec3));

  % === Extract and process output channel data ===

  cd(OutputFilesDir);
  nch = 1;
  fid = 1;
  while ~eq(fid,-1)
      ch_filename = sprintf('Output_Channels_%d.dat',nch);

      fid = fopen(ch_filename,'r');
      if ~eq(fid,-1) % i.e. channel file exists
          ldata = [];
          vdata = [];
          qxdata = [];
          qydata = [];
          qzdata = [];
          bdata = [];
          vxdata = [];
          vydata = [];
          vzdata = [];
          wdata = [];
          nchdata = [];

          fprintf(1,'++ Channel data file name = %s \n',ch_filename);

          tmp = importdata(ch_filename,' ',1);
          if isfield(tmp,'data') % i.e. channel data exist
              % fprintf(1,'%4d  %4d  %4d  \n', nch, size(tmp.data,1), fid);
              [status, result] = system(sprintf('head -1 %s | awk ''{print $5}''',ch_filename));
              wch = str2num(result); % channel frequency    
              ldata = [ldata; tmp.data(:,1)];
              vdata = [vdata; tmp.data(:,2)];
              qxdata = [qxdata; tmp.data(:,3)];
              qydata = [qydata; tmp.data(:,4)];
              qzdata = [qzdata; tmp.data(:,5)];
              bdata = [bdata; tmp.data(:,6)];
              vxdata = [vxdata; tmp.data(:,7)];
              vydata = [vydata; tmp.data(:,8)];
              vzdata = [vzdata; tmp.data(:,9)];
              wdata = [wdata; wch*ones(size(tmp.data,1),1)];
              nchdata = [nchdata; nch*ones(size(tmp.data,1),1)];

              fprintf(1,'\t No. of scattering (incoming+outgoing) channels = %d \n',numel(wdata));


              % === Reprocess output channel data to find correct band indices ===
              tic;
            
              for nq = 1:numel(ldata)
                  qvec = [qxdata(nq); qydata(nq); qzdata(nq)];
                  vqvec = [vxdata(nq); vydata(nq); vzdata(nq)];  
            
                  if ldata(nq)<0
                      [wwq,Uq,Hq,vq] = GetBranchSortedFreq(LeftPhon.S,qvec,LeftPhon.qcutoff_G);    
                  else
                      [wwq,Uq,Hq,vq] = GetBranchSortedFreq(RightPhon.S,qvec,RightPhon.qcutoff_G);    
                  end
    
                  wwq = diag(wwq);
                  [y_tmp,ind_tmp] = sort(abs(wdata(nq)^2-wwq));
                  wdegeneracy = numel(wwq(abs((wdata(nq)-sqrt(wwq))/wdata(nq)) <delta_w)); % no. of modes at this q-point with same freq as mode nq

                  if wdegeneracy>1 % i.e. the mode is degenerate
                      nrange = (ldata(1:nq,1)==ldata(nq,1))&(vdata(1:nq,1)==vdata(nq,1))&(wdata(1:nq,1)==wdata(nq,1)); % same lead/vel sign/frequency, 
                      subdqxdata = qxdata(nrange) - qxdata(nq);
                      subdqydata = qydata(nrange) - qydata(nq);
                      subdqzdata = qzdata(nrange) - qzdata(nq);
                      error_qdata = sqrt(subdqxdata.^2 + subdqydata.^2 + subdqzdata.^2)/sqrt(qxdata(nq)^2 + qydata(nq)^2 + qzdata(nq)^2);
                      subdvxdata = vxdata(nrange) - vxdata(nq);
                      subdvydata = vydata(nrange) - vydata(nq);
                      subdvzdata = vzdata(nrange) - vzdata(nq);
                      error_vdata = sqrt(subdvxdata.^2 + subdvydata.^2 + subdvzdata.^2)/sqrt(vxdata(nq)^2 + vydata(nq)^2 + vzdata(nq)^2);
                      ndegen = numel(error_qdata((error_qdata<delta_q)&(error_vdata<delta_v))); % no. of previous incoming modes with same q and vel at this freq
                      % ndegen = numel(error_qdata(error_qdata<delta_q)); % no. of previous incoming modes with same q at this freq

                      for n_tmp = 1:numel(ind_tmp)
                          error_w = abs(wdata(nq)-sqrt(wwq(ind_tmp(n_tmp),1)))/wdata(nq);
                          error_v = norm(vqvec-vq(:,ind_tmp(n_tmp)))/norm(vqvec);
                          if (error_w<delta_w)&(error_v<delta_v)
                              ndegen = ndegen - 1;
                          end
            
                          if ndegen==0
                              bdata(nq,1) = ind_tmp(n_tmp);
                              break;
                          end
                      end
                  else
                      bdata(nq,1) = ind_tmp(1);
                  end
                
                  if or(mod(nq,200)==0,nq==numel(ldata))
                      fprintf(1,'\t No. of channels processed = %6d, time elapsed = %12.4f \n', nq, toc);
                  end
              end
          end

          fclose(fid);

          % === Write reprocessed output channel data to file ===

          [status, header_text] = system(sprintf('head -1 %s | awk ''{print $0}''',ch_filename));
          fid = fopen(ch_filename,'w');         

          fprintf(fid,'%s',header_text);
          for n = 1:numel(wdata)
              if nchdata(n,1)==nch
                  fprintf(fid,'%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e \n', ...
                      ldata(n,1), vdata(n,1), qxdata(n,1), qydata(n,1), qzdata(n,1), ...
                      bdata(n,1), vxdata(n,1), vydata(n,1), vzdata(n,1));
              end
          end
          if and(fid~=1,fid~=-1)
              fclose(fid);   
          end

          nch = nch + 1;
      end
  end

  nchmax = nch - 1;
  fprintf(1,'++ No. of channel data files processed = %d \n',nchmax);

  cd(CurrDir);

  % === Extract transmission data ===

  wdata_tr = [];
  bdata_tr = [];
  qxdata_tr = [];
  qydata_tr = [];
  qzdata_tr = [];
  cxdata_tr = [];
  pdata_tr = [];
  tdata_tr = [];
  vxdata_tr = [];
  vydata_tr = [];
  vzdata_tr = [];
  ldata_tr = [];

  clear tmp;

  cd(OutputFilesDir);
  tr_filename_left = sprintf('Output_Left_Transmission.dat');
  tmp.data = load(tr_filename_left);
  if isfield(tmp,'data') % i.e. data exist
      wdata_tr  = tmp.data(:,1); 
      bdata_tr  = tmp.data(:,2); 
      qxdata_tr = tmp.data(:,3); 
      qydata_tr = tmp.data(:,4); 
      qzdata_tr = tmp.data(:,5); 
      cxdata_tr = tmp.data(:,6); 
      pdata_tr  = tmp.data(:,7); 
      tdata_tr  = tmp.data(:,8);
      vxdata_tr = tmp.data(:,9);
      vydata_tr = tmp.data(:,10);
      vzdata_tr = tmp.data(:,11);
      ldata_tr  = -1*ones(size(tmp.data(:,1)));
  end

  tr_filename_right = sprintf('Output_Right_Transmission.dat');
  tmp.data = load(tr_filename_right);
  if isfield(tmp,'data') % i.e. data exist
      wdata_tr  = [wdata_tr; tmp.data(:,1)]; 
      bdata_tr  = [bdata_tr; tmp.data(:,2)]; 
      qxdata_tr = [qxdata_tr; tmp.data(:,3)]; 
      qydata_tr = [qydata_tr; tmp.data(:,4)]; 
      qzdata_tr = [qzdata_tr; tmp.data(:,5)]; 
      cxdata_tr = [cxdata_tr; tmp.data(:,6)]; 
      pdata_tr  = [pdata_tr; tmp.data(:,7)]; 
      tdata_tr  = [tdata_tr; tmp.data(:,8)];
      vxdata_tr = [vxdata_tr; tmp.data(:,9)];
      vydata_tr = [vydata_tr; tmp.data(:,10)];
      vzdata_tr = [vzdata_tr; tmp.data(:,11)];
      ldata_tr  = [ldata_tr; +1*ones(size(tmp.data(:,1)))];
  end
  cd(CurrDir);

  fprintf(1,'++ No. of left lead transmission channels = %d \n',numel(bdata_tr(ldata_tr<0)));
  fprintf(1,'++ No. of right lead transmission channels = %d \n',numel(bdata_tr(ldata_tr>0)));

  % === Reprocess transmission data to find correct band indices ===
  tic;

  for nq = 1:numel(bdata_tr)
      qvec = [qxdata_tr(nq); qydata_tr(nq); qzdata_tr(nq)]; 
      vqvec = [vxdata_tr(nq); vydata_tr(nq); vzdata_tr(nq)]; 

      if ldata_tr(nq)<0
          [wwq,Uq,Hq,vq] = GetBranchSortedFreq(LeftPhon.S,qvec,LeftPhon.qcutoff_G);    
      else
          [wwq,Uq,Hq,vq] = GetBranchSortedFreq(RightPhon.S,qvec,RightPhon.qcutoff_G);    
      end
    
      wwq = diag(wwq);
      [y_tmp,ind_tmp] = sort(abs(wdata_tr(nq)^2-wwq));

      wdegeneracy = numel(wwq(abs((wdata_tr(nq)-sqrt(wwq))/wdata_tr(nq))<delta_w)); % no. of modes at this q-point with same freq as mode nq

      if wdegeneracy>1 % i.e. the mode is degenerate
          nrange = (ldata_tr(1:nq,1)==ldata_tr(nq,1))&(wdata_tr(1:nq,1)==wdata_tr(nq,1));
          subdqxdata_tr = qxdata_tr(nrange) - qxdata_tr(nq);
          subdqydata_tr = qydata_tr(nrange) - qydata_tr(nq);
          subdqzdata_tr = qzdata_tr(nrange) - qzdata_tr(nq);
          error_qdata_tr = sqrt(subdqxdata_tr.^2 + subdqydata_tr.^2 + subdqzdata_tr.^2)/sqrt(qxdata_tr(nq)^2 + qydata_tr(nq)^2 + qzdata_tr(nq)^2);
          subdvxdata_tr = vxdata_tr(nrange) - vxdata_tr(nq);
          subdvydata_tr = vydata_tr(nrange) - vydata_tr(nq);
          subdvzdata_tr = vzdata_tr(nrange) - vzdata_tr(nq);
          error_vdata_tr = sqrt(subdvxdata_tr.^2 + subdvydata_tr.^2 + subdvzdata_tr.^2)/sqrt(vxdata_tr(nq)^2 + vydata_tr(nq)^2 + vzdata_tr(nq)^2);
          ndegen = numel(error_qdata_tr((error_qdata_tr<delta_q)&(error_vdata_tr<delta_v))); % no. of previous incoming modes with same q and vel at this freq
          % ndegen = numel(error_qdata_tr(error_qdata_tr<delta_q)); % no. of previous incoming modes with same q at this freq

          for n_tmp = 1:numel(ind_tmp)
              error_w = abs(wdata_tr(nq)-sqrt(wwq(ind_tmp(n_tmp),1)))/wdata_tr(nq);
              error_v = norm(vqvec-vq(:,ind_tmp(n_tmp)))/norm(vqvec);
              if (error_w<delta_w)&(error_v<delta_v)
                  ndegen = ndegen - 1;
              end

              if ndegen==0
                  bdata_tr(nq,1) = ind_tmp(n_tmp);
                  break;
              end
          end

          %{ 
          fprintf(1,'\n'); % DEBUG
          fprintf(1,'\t %3d : %14.6e : %14.6e %14.6e %14.6e : %3d %3d \n',nq, wdata_tr(nq), qvec(1), qvec(2), qvec(3), ind_tmp(1), bdata_tr(nq,1)); % DEBUG
          % fprintf(1,'\t %14.6e',sqrt(wwq')); % DEBUG
          fprintf(1,'\n'); % DEBUG
          fprintf(1,'\t\t %14.6e %14.6e %14.6e \n',vxdata_tr(nq),vydata_tr(nq),vzdata_tr(nq)); % DEBUG
          fprintf(1,'\t\t %14.6e %14.6e %14.6e \n',vq(1,bdata_tr(nq,1)),vq(2,bdata_tr(nq,1)),vq(3,bdata_tr(nq,1))); % DEBUG
          %}
      else
          bdata_tr(nq,1) = ind_tmp(1);
      end
    
      if or(mod(nq,200)==0,nq==numel(ldata_tr))
          fprintf(1,'\t No. of transmission channels processed = %6d, time elapsed = %12.4f \n', nq, toc);
      end
  end
  % return; % DEBUG	

  % === Write reprocessed transmission data ===
  cd(OutputFilesDir);
  fid_left = fopen(tr_filename_left,'w');
  fid_right = fopen(tr_filename_right,'w');
  for nq = 1:numel(bdata_tr)
      if ldata_tr(nq)==-1
          fprintf(fid_left,'%12.4e %3d %12.4e %12.4e %12.4e %12.4e %14.6e %14.6e %12.4e %12.4e %12.4e \n', ...
              wdata_tr(nq), bdata_tr(nq), qxdata_tr(nq), qydata_tr(nq), qzdata_tr(nq), cxdata_tr(nq), ... 
              pdata_tr(nq), tdata_tr(nq), vxdata_tr(nq), vydata_tr(nq), vzdata_tr(nq));            
      elseif ldata_tr(nq)==1
          fprintf(fid_right,'%12.4e %3d %12.4e %12.4e %12.4e %12.4e %14.6e %14.6e %12.4e %12.4e %12.4e \n', ...
              wdata_tr(nq), bdata_tr(nq), qxdata_tr(nq), qydata_tr(nq), qzdata_tr(nq), cxdata_tr(nq), ... 
              pdata_tr(nq), tdata_tr(nq), vxdata_tr(nq), vydata_tr(nq), vzdata_tr(nq));            
      end
  end
  fclose(fid_left);
  fclose(fid_right);
  cd(CurrDir);
  %{

  %}
end

% =========================================================================
% =========================================================================

function [wwq_out,Uq_out,Hq_out,vq_out] = GetBranchSortedFreq(S_in,qvec_in,qcutoff_in)
  if norm(qvec_in)<=qcutoff_in
      nstepmax = 1;
  else
      nstepmax = 200;
  end

  for nstep = 1:nstepmax
      if norm(qvec_in)<=qcutoff_in
          qvec = qvec_in;
      else
          qvec_i = qvec_in*qcutoff_in/norm(qvec_in);
          qvec_f = qvec_in;
          qvec = qvec_i + nstep*1.0/nstepmax*(qvec_f-qvec_i);
      end
    
      Hq = zeros(size(S_in(1).H));
      for n = 1:numel(S_in)
          Hq = Hq + S_in(n).H * exp(1i*S_in(n).rvec'*qvec);
      end
      Hq = 0.5*(Hq + Hq');
      [Uq, wwq] = eig(Hq);
    
      if nstep==1
          [wwq_tmp,ind_tmp] = sort(diag(real(wwq)));
          wwq = diag(wwq_tmp); % sorted eigenfreq squared
          Uq = Uq(:,ind_tmp); % sorted eigenvec
          
          Uq_old = Uq;
          Hq_old = Hq;
          wwq_old = wwq;
      else
          Uq_new = Uq;
          V = Hq - Hq_old;
          wwq0 = wwq_old;
          
          for n = 1:size(wwq0,1)
              dUq(:,n) = Uq_old*diag(pinv(diag(wwq0(n) - wwq0)))*Uq_old'*V*Uq_old(:,n);
              Uq1_old(:,n) = Uq_old(:,n) + dUq(:,n); % 'shifted' eigenvec from perturbation theory
              Uq1_old(:,n) = Uq1_old(:,n)/norm(Uq1_old(:,n));
          end
          Sq = abs(Uq1_old'*Uq_new);
          
          [y_tmp, ind_tmp] = max(Sq,[],2);
          Uq = Uq(:,ind_tmp); % exact eigenvec reordered to previous eigenvec
          wwq_tmp = real(diag(wwq));
          wwq = diag(wwq_tmp(ind_tmp)); % exact eigenfreq squared reordered to previous eigenvec
          
          Uq_old = Uq;
          Hq_old = Hq;
          wwq_old = wwq;
      end
  end

  wwq_out = wwq; % ordered eigenfreq squared (as diagonal matrix)
  Uq_out = Uq;   % ordered eigenvec
  Hq_out = Hq;   % 'Hamiltonian'

  %{
  %}
  Hq_x = S_in(1).H*exp(1i*S_in(1).rvec'*qvec_in)*1i*S_in(1).rvec(1,1);
  Hq_y = S_in(1).H*exp(1i*S_in(1).rvec'*qvec_in)*1i*S_in(1).rvec(2,1);
  Hq_z = S_in(1).H*exp(1i*S_in(1).rvec'*qvec_in)*1i*S_in(1).rvec(3,1);
  for np = 2:numel(S_in)
      Hq_x = Hq_x + S_in(np).H*exp(1i*S_in(np).rvec'*qvec_in) ...
             *1i*S_in(np).rvec(1,1);
      Hq_y = Hq_y + S_in(np).H*exp(1i*S_in(np).rvec'*qvec_in) ...
             *1i*S_in(np).rvec(2,1);
      Hq_z = Hq_z + S_in(np).H*exp(1i*S_in(np).rvec'*qvec_in) ...
             *1i*S_in(np).rvec(3,1);
  end

  vq_out = zeros(3,size(wwq,1));

  for nq = 1:size(wwq,1)
      wq = sqrt(wwq(nq,nq));
      vq_out(1,nq) = real(1.0/(2.0*wq)*Uq(:,nq)'*Hq_x*Uq(:,nq));
      vq_out(2,nq) = real(1.0/(2.0*wq)*Uq(:,nq)'*Hq_y*Uq(:,nq));
      vq_out(3,nq) = real(1.0/(2.0*wq)*Uq(:,nq)'*Hq_z*Uq(:,nq));
  end
end



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

