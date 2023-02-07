function WritePhononScatteringOutput(DataFilesDir,wvec_scatter,wvec,LeftPhonon,RightPhonon,PhononData)
% NOTE: This function writes:
% NOTE:   (1) the transpose of the transition probability matrix or W^T ("Out_Scattering_%d.dat")
% NOTE:   (2) the Hermitean transpose of the S matrix or S^dagger ("Output_SMatrix_%d.dat")
% NOTE:   (3) the phonon channel data describing the order of the row and column indices ("Output_Channels_%d.dat")
% NOTE:   (4) the specularity data ("Output_Specularity_%d.dat")
% NOTE: Format of output files:
% NOTE:   (1) NxN matrix [W^T], with each row (column) corresponding to an incoming (outgoing) phonon channel 
% NOTE:   (2) Nx2N matrix [Real(S^dagger) Imag(S^dagger)], with each row (column) corresponding to an incoming (outgoing) phonon channel 
% NOTE:   (3) 9 columns: Lead, velocity direction, k_x, k_y, k_z, band index, v_x, v_y, v_z 
% NOTE:   (4) 10 columns: Lead, velocity direction, k_x, k_y, k_z, band index, v_x, v_y, v_z, specularity probability 

  CurrDir = pwd; % current directory containing AGF code .m files

  cd(DataFilesDir); % change to AGF output directory

  if eq(wvec_scatter(1,1),0)
      warning('No scattering data to be written.');
      exit;
  else
      % disp(numel(wvec_scatter));
      nwrange = 1:numel(wvec);

      for nw = 1:numel(wvec_scatter)
          w_scatter = wvec_scatter(nw); % scattering frequency
          nw_scatter = nwrange(eq(wvec,w_scatter)); % index position of scattering freq. in list of input freqs.

          if gt(numel(nw_scatter),0)
              % === LEFT-TO-RIGHT ===
              b_in_L_temp  = [LeftPhonon(nw_scatter).VecB_minus_adv]; % band index of incident phonons in left lead
              b_out_R_temp = [RightPhonon(nw_scatter).VecB_plus];     % band index of outgoing phonons in right lead
              b_out_L_temp = [LeftPhonon(nw_scatter).VecB_minus];     % band index of outgoing phonons in left lead

              nrange_in_L = gt(b_in_L_temp,0);   
              nrange_out_R = gt(b_out_R_temp,0);
              nrange_out_L = gt(b_out_L_temp,0);

              tmatrix_tmp = PhononData(nw_scatter).t_RL;
              rmatrix_tmp = PhononData(nw_scatter).r_LL;
              tmatrix_tmp = tmatrix_tmp(nrange_out_R,nrange_in_L);
              rmatrix_tmp = rmatrix_tmp(nrange_out_L,nrange_in_L);
              Smatrix_L = [rmatrix_tmp' tmatrix_tmp'];
              % tmatrix_L = [tmatrix_tmp; rmatrix_tmp]';
              % tmatrix = [tmatrix sum(abs(tmatrix).^2,2)];

              % === RIGHT-TO-LEFT ===
              b_in_R_temp  = [RightPhonon(nw_scatter).VecB_plus_adv]; % band index of incident phonons in right lead
              b_out_L_temp = [LeftPhonon(nw_scatter).VecB_minus];     % band index of outgoing phonons in left lead
              b_out_R_temp = [RightPhonon(nw_scatter).VecB_plus];     % band index of outgoing phonons in right lead

              nrange_in_R = gt(b_in_R_temp,0);   
              nrange_out_L = gt(b_out_L_temp,0);
              nrange_out_R = gt(b_out_R_temp,0);

              tmatrix_tmp = PhononData(nw_scatter).t_LR;
              rmatrix_tmp = PhononData(nw_scatter).r_RR;
              tmatrix_tmp = tmatrix_tmp(nrange_out_L,nrange_in_R);
              rmatrix_tmp = rmatrix_tmp(nrange_out_R,nrange_in_R);
              Smatrix_R = [tmatrix_tmp' rmatrix_tmp'];
              % tmatrix_R = [rmatrix_tmp; tmatrix_tmp]';

              % === write scattering (transition probability) data ===
              filename_scatter = sprintf('Output_Scattering_%d.dat',nw);
              fprintf(1,'\t  <%s> \n', filename_scatter);
              fid = fopen(filename_scatter,'w');
              fprintf(fid, '# Angular frequency = %16.8e \n', w_scatter);

              for n = 1:size(Smatrix_L,1) 
                  fprintf(fid, '%14.6e', abs(Smatrix_L(n,:)).^2);
                  fprintf(fid, '\n');
              end

              for n = 1:size(Smatrix_R,1) 
                  fprintf(fid, '%14.6e', abs(Smatrix_R(n,:)).^2);
                  fprintf(fid, '\n');
              end

              fclose(fid);


              % === write S-matrix data (NEW!!!) ===
              filename_smatrix = sprintf('Output_SMatrix_%d.dat',nw);
              fprintf(1,'\t  <%s> \n', filename_smatrix);
              fid = fopen(filename_smatrix,'w');
              fprintf(fid, '# Angular frequency = %16.8e \n', w_scatter);

              for n = 1:size(Smatrix_L,1) 
                  fprintf(fid, '%14.6e', real(Smatrix_L(n,:)));
                  fprintf(fid, '%14.6e', imag(Smatrix_L(n,:)));
                  fprintf(fid, '\n');
              end

              for n = 1:size(Smatrix_R,1) 
                  fprintf(fid, '%14.6e', real(Smatrix_R(n,:)));
                  fprintf(fid, '%14.6e', imag(Smatrix_R(n,:)));
                  fprintf(fid, '\n');
              end

              fclose(fid);


              % --- wave vector and group velocity data of incident phonons in left lead ---
              qx_in_L_temp = [LeftPhonon(nw_scatter).VecQx_minus_adv]; % wave vector (x-dir) of incident phonons  
              qy_in_L_temp = [LeftPhonon(nw_scatter).VecQy_minus_adv]; % wave vector (y-dir) of incident phonons  
              qz_in_L_temp = [LeftPhonon(nw_scatter).VecQz_minus_adv]; % wave vector (z-dir) of incident phonons  
              velx_in_L_temp = [LeftPhonon(nw_scatter).VecVelx_minus_adv]; % group velocity (x-dir) of incident phonons  
              vely_in_L_temp = [LeftPhonon(nw_scatter).VecVely_minus_adv]; % group velocity (y-dir) of incident phonons  
              velz_in_L_temp = [LeftPhonon(nw_scatter).VecVelz_minus_adv]; % group velocity (z-dir) of incident phonons  
    
              qx_in_L_temp = qx_in_L_temp(nrange_in_L);
              qy_in_L_temp = qy_in_L_temp(nrange_in_L);
              qz_in_L_temp = qz_in_L_temp(nrange_in_L);        
              velx_in_L_temp = velx_in_L_temp(nrange_in_L);
              vely_in_L_temp = vely_in_L_temp(nrange_in_L);
              velz_in_L_temp = velz_in_L_temp(nrange_in_L);
              b_in_L_temp  = b_in_L_temp(nrange_in_L);

              % --- wave vector and group velocity data of incident phonons in right lead ---
              qx_in_R_temp = [RightPhonon(nw_scatter).VecQx_plus_adv]; % wave vector (x-dir) of incident phonons  
              qy_in_R_temp = [RightPhonon(nw_scatter).VecQy_plus_adv]; % wave vector (y-dir) of incident phonons  
              qz_in_R_temp = [RightPhonon(nw_scatter).VecQz_plus_adv]; % wave vector (z-dir) of incident phonons  
              velx_in_R_temp = [RightPhonon(nw_scatter).VecVelx_plus_adv]; % group velocity (x-dir) of incident phonons  
              vely_in_R_temp = [RightPhonon(nw_scatter).VecVely_plus_adv]; % group velocity (y-dir) of incident phonons  
              velz_in_R_temp = [RightPhonon(nw_scatter).VecVelz_plus_adv]; % group velocity (z-dir) of incident phonons  
    
              qx_in_R_temp = qx_in_R_temp(nrange_in_R);
              qy_in_R_temp = qy_in_R_temp(nrange_in_R);
              qz_in_R_temp = qz_in_R_temp(nrange_in_R);        
              velx_in_R_temp = velx_in_R_temp(nrange_in_R);
              vely_in_R_temp = vely_in_R_temp(nrange_in_R);
              velz_in_R_temp = velz_in_R_temp(nrange_in_R);
              b_in_R_temp  = b_in_R_temp(nrange_in_R);

              % --- wave vector and group velocity data of outgoing phonons in right lead ---
              qx_out_R_temp = [RightPhonon(nw_scatter).VecQx_plus]; % wave vector (x-dir) of outgoing phonons
              qy_out_R_temp = [RightPhonon(nw_scatter).VecQy_plus]; % wave vector (y-dir) of outgoing phonons
              qz_out_R_temp = [RightPhonon(nw_scatter).VecQz_plus]; % wave vector (z-dir) of outgoing phonons    
              velx_out_R_temp = [RightPhonon(nw_scatter).VecVelx_plus]; % group velocity (x-dir) of outgoing phonons
              vely_out_R_temp = [RightPhonon(nw_scatter).VecVely_plus]; % group velocity (y-dir) of outgoing phonons
              velz_out_R_temp = [RightPhonon(nw_scatter).VecVelz_plus]; % group velocity (z-dir) of outgoing phonons

              qx_out_R_temp = qx_out_R_temp(nrange_out_R);
              qy_out_R_temp = qy_out_R_temp(nrange_out_R);
              qz_out_R_temp = qz_out_R_temp(nrange_out_R);
              velx_out_R_temp = velx_out_R_temp(nrange_out_R);
              vely_out_R_temp = vely_out_R_temp(nrange_out_R);
              velz_out_R_temp = velz_out_R_temp(nrange_out_R);
              b_out_R_temp  = b_out_R_temp(nrange_out_R);

              % --- wave vector and group velocity data of outgoing phonons in left lead ---
              qx_out_L_temp = [LeftPhonon(nw_scatter).VecQx_minus]; % wave vector (x-dir) of outgoing phonons  
              qy_out_L_temp = [LeftPhonon(nw_scatter).VecQy_minus]; % wave vector (y-dir) of outgoing phonons  
              qz_out_L_temp = [LeftPhonon(nw_scatter).VecQz_minus]; % wave vector (z-dir) of outgoing phonons  
              velx_out_L_temp = [LeftPhonon(nw_scatter).VecVelx_minus]; % group (x-dir) of outgoing phonons  
              vely_out_L_temp = [LeftPhonon(nw_scatter).VecVely_minus]; % group (y-dir) of outgoing phonons  
              velz_out_L_temp = [LeftPhonon(nw_scatter).VecVelz_minus]; % group (z-dir) of outgoing phonons  

              qx_out_L_temp = qx_out_L_temp(nrange_out_L);
              qy_out_L_temp = qy_out_L_temp(nrange_out_L);
              qz_out_L_temp = qz_out_L_temp(nrange_out_L);
              velx_out_L_temp = velx_out_L_temp(nrange_out_L);
              vely_out_L_temp = vely_out_L_temp(nrange_out_L);
              velz_out_L_temp = velz_out_L_temp(nrange_out_L);
              b_out_L_temp  = b_out_L_temp(nrange_out_L);

              % ===== (rationalize wave vectors by 2 pi) =====
              qx_in_L_temp = 2*pi*qx_in_L_temp;
              qy_in_L_temp = 2*pi*qy_in_L_temp;
              qz_in_L_temp = 2*pi*qz_in_L_temp;
              qx_in_R_temp = 2*pi*qx_in_R_temp;
              qy_in_R_temp = 2*pi*qy_in_R_temp;
              qz_in_R_temp = 2*pi*qz_in_R_temp;
        
              qx_out_L_temp = 2*pi*qx_out_L_temp;
              qy_out_L_temp = 2*pi*qy_out_L_temp;
              qz_out_L_temp = 2*pi*qz_out_L_temp;
              qx_out_R_temp = 2*pi*qx_out_R_temp;
              qy_out_R_temp = 2*pi*qy_out_R_temp;
              qz_out_R_temp = 2*pi*qz_out_R_temp;

              % === write channel data === 
              filename_channel = sprintf('Output_Channels_%d.dat',nw);
              fprintf(1,'\t  <%s> \n', filename_channel);
              fid = fopen(filename_channel,'w');
              fprintf(fid, '# Angular frequency = %16.8e : (lead, v_x sign, k_x, k_y, k_z, band index, v_x, v_y, v_z ) \n', w_scatter);

              % --- incident phonons in left lead (type -1,+1) ---
              for n = 1:numel(b_in_L_temp) 
                  fprintf(fid, '%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e \n', ...
                    -1, +1, qx_in_L_temp(n), qy_in_L_temp(n), qz_in_L_temp(n), ...
                    b_in_L_temp(n), velx_in_L_temp(n), vely_in_L_temp(n), velz_in_L_temp(n));
              end

              % --- incident phonons in right lead (type -1,+1) ---
              for n = 1:numel(b_in_R_temp) 
                  fprintf(fid, '%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e \n', ...
                    +1, -1, qx_in_R_temp(n), qy_in_R_temp(n), qz_in_R_temp(n), ...
                    b_in_R_temp(n), velx_in_R_temp(n), vely_in_R_temp(n), velz_in_R_temp(n));
              end

              % --- outgoing phonons in left lead (type -1,-1) ---
              for n = 1:numel(b_out_L_temp) 
                  fprintf(fid, '%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e \n', ...
                    -1, -1, qx_out_L_temp(n), qy_out_L_temp(n), qz_out_L_temp(n), ...
                    b_out_L_temp(n), velx_out_L_temp(n), vely_out_L_temp(n), velz_out_L_temp(n));
              end

              % --- outgoing phonons in right lead (type +1,+1) ---
              for n = 1:numel(b_out_R_temp) 
                  fprintf(fid, '%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e \n', ...
                    +1, +1, qx_out_R_temp(n), qy_out_R_temp(n), qz_out_R_temp(n), ...
                    b_out_R_temp(n), velx_out_R_temp(n), vely_out_R_temp(n), velz_out_R_temp(n));
              end
              fclose(fid);


              % === write channel phonon specularity parameter data === 
              filename_specular = sprintf('Output_Specularity_%d.dat',nw);
              fprintf(1,'\t  <%s> \n', filename_specular);
              fid = fopen(filename_specular,'w');
              fprintf(fid, '# Angular frequency = %16.8e \n', w_scatter);
              
              % --- incident phonons in left lead (type -1,+1) ---
              rmatrix_tmp = PhononData(nw_scatter).r_LL;
              rmatrix_tmp = rmatrix_tmp(nrange_out_L,nrange_in_L);

              for n1 = 1:numel(b_in_L_temp) 
                  qx1 = qx_in_L_temp(n1);
                  qy1 = qy_in_L_temp(n1);
                  qz1 = qz_in_L_temp(n1);
                  b1 = b_in_L_temp(n1);
                  vx1 = velx_in_L_temp(n1);
                  vy1 = vely_in_L_temp(n1);
                  vz1 = velz_in_L_temp(n1);
                  qnorm1 = norm([qx1 qy1 qz1]);
                 
                  for n2 = 1:numel(b_out_L_temp) 
                      qx2 = qx_out_L_temp(n2);
                      qy2 = qy_out_L_temp(n2);
                      qz2 = qz_out_L_temp(n2);
                      % qnorm2 = norm([qx1 qy1 qz1]);
                      delta_qnorm = norm([(qx1+qx2) (qy1-qy2) (qz1-qz2)]);
                      if lt(delta_qnorm/qnorm1,1E-3)
                          p1_spec = abs(rmatrix_tmp(n2,n1)).^2;
                      end
                  end
                  fprintf(fid, '%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e %14.6e \n', ...
                    -1, +1, qx1, qy1, qz1, b1, vx1, vy1, vz1, p1_spec);
              end

              rmatrix_tmp = PhononData(nw_scatter).r_RR;
              rmatrix_tmp = rmatrix_tmp(nrange_out_R,nrange_in_R);

              for n1 = 1:numel(b_in_R_temp) 
                  qx1 = qx_in_R_temp(n1);
                  qy1 = qy_in_R_temp(n1);
                  qz1 = qz_in_R_temp(n1);
                  b1 = b_in_R_temp(n1);
                  vx1 = velx_in_R_temp(n1);
                  vy1 = vely_in_R_temp(n1);
                  vz1 = velz_in_R_temp(n1);
                  qnorm1 = norm([qx1 qy1 qz1]);
                 
                  for n2 = 1:numel(b_out_R_temp) 
                      qx2 = qx_out_R_temp(n2);
                      qy2 = qy_out_R_temp(n2);
                      qz2 = qz_out_R_temp(n2);
                      % qnorm2 = norm([qx1 qy1 qz1]);
                      delta_qnorm = norm([(qx1+qx2) (qy1-qy2) (qz1-qz2)]);
                      if lt(delta_qnorm/qnorm1,1E-3)
                          p1_spec = abs(rmatrix_tmp(n2,n1)).^2;
                      end
                  end
                  fprintf(fid, '%3d %3d %14.6e %14.6e %14.6e %3d %14.6e %14.6e %14.6e %14.6e \n', ...
                    +1, -1, qx1, qy1, qz1, b1, vx1, vy1, vz1, p1_spec);
              end

              fclose(fid); 

          end
      end
  end

  cd(CurrDir); % return to directory containing AGF code .m files

end
