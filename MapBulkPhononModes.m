function LeadPhonon = MapBulkPhononModes(w_in,LeadPhonon_in,Lead,LeadPhonDisp)
% NOTE: This function maps the AGF modes to their bulk lattice phonon modes 
% NOTE: and determines their primitive BZ wave vector and group velocity
% NOTE: w_in = input frequency
% NOTE: LeadPhonon_in = data structure storing AGF lead phonon eigenvectors, velocities, wave vectors, etc
% NOTE: Lead = data structure storing force constant matrices and transverse cell info for lead
% NOTE: LeadPhonDisp = data structure for bulk lattice phonon dispersion
% NOTES: This functions uses "GenerateImageBZPoints.m"

  w = w_in;
  LeadPhonon = LeadPhonon_in;

  % return; % DEBUG
  epsilon_cutoff_w = 1E-4; % matching cutoff parameter for frequency
  epsilon_cutoff_q = 1E-4; % matching cutoff parameter for wave vector
  % nwmax = numel(wvec);

  % if and(eq(Lead.n_tran1,1),eq(Lead.n_tran2,1))
  %   fprintf(1,'<DEBUG> System is quasi-1D \n');
  % end

  % ------------------------------------------------------------------------------
  % Find unfolded wave vectors 
  %   Axis 1 = longitudinal direction
  %   Axis 2 = transverse direction 1
  %   Axis 3 = transverse direction 2
  % ------------------------------------------------------------------------------

  Gvec1 = Lead.gvec_long;  % reciprocal lattice vector (AGF, non-primitive)
  Gvec2 = Lead.gvec_tran1; % reciprocal lattice vector (AGF, non-primitive)
  Gvec3 = Lead.gvec_tran2; % reciprocal lattice vector (AGF, non-primitive)
  PrimGvec1 = LeadPhonDisp.Gvec1; % reciprocal lattice vector (Brillouin Zone, primitive)
  PrimGvec2 = LeadPhonDisp.Gvec2; % reciprocal lattice vector (Brillouin Zone, primitive)
  PrimGvec3 = LeadPhonDisp.Gvec3; % reciprocal lattice vector (Brillouin Zone, primitive)
  
  PrimCell = LeadPhonDisp.PrimCell; % contains parameters for computing primitive phonon dispersion 
  
  ScaleFactor = abs(PrimGvec1*cross(PrimGvec2,PrimGvec3)')/abs(Gvec1*cross(Gvec2,Gvec3)');
    % ratio of primitive BZ volume to AGF BZ volume, should be an integer >= 1; also equal number of image points 

  % =================================================================================================
  % Retarded plus phonon channels
  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;   % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;   % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_plus]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_plus]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;  % set value between -0.5 and 0.4999...
  p_temp = gt(abs(qfrac1_temp),1E-5); % indices of propagating phonons

  U_temp = LeadPhonon.U_plus;
  [qdata_temp, Udata_temp, vdata_temp, bdata_temp] = ...
    GenerateImageBZPoints(w,vlong_temp,qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,U_temp,Lead,LeadPhonDisp); % DEBUG

  LeadPhonon.U_plus = Udata_temp;
  LeadPhonon.VecQx_plus = qdata_temp(:,1);
  LeadPhonon.VecQy_plus = qdata_temp(:,2);
  LeadPhonon.VecQz_plus = qdata_temp(:,3);
  LeadPhonon.VecVelx_plus = vdata_temp(:,1);
  LeadPhonon.VecVely_plus = vdata_temp(:,2);
  LeadPhonon.VecVelz_plus = vdata_temp(:,3);
  LeadPhonon.VecB_plus  = bdata_temp(:,1);

  % =================================================================================================
  % Advanced plus phonon channels
  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_plus_adv]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_plus_adv]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;           % set value between -0.5 and 0.4999...
  p_temp = gt(abs(qfrac1_temp),1E-5); % indices of propagating phonons

  U_temp = LeadPhonon.U_plus_adv;
  [qdata_temp, Udata_temp, vdata_temp, bdata_temp] = ...
    GenerateImageBZPoints(w,vlong_temp,qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,U_temp,Lead,LeadPhonDisp); % DEBUG

  LeadPhonon.U_plus_adv = Udata_temp;
  LeadPhonon.VecQx_plus_adv = qdata_temp(:,1);
  LeadPhonon.VecQy_plus_adv = qdata_temp(:,2);
  LeadPhonon.VecQz_plus_adv = qdata_temp(:,3);
  LeadPhonon.VecVelx_plus_adv = vdata_temp(:,1);
  LeadPhonon.VecVely_plus_adv = vdata_temp(:,2);
  LeadPhonon.VecVelz_plus_adv = vdata_temp(:,3);
  LeadPhonon.VecB_plus_adv  = bdata_temp(:,1);

  % =================================================================================================
  % Retarded minus phonon channels
  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_minus]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_minus]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;           % set value between -0.5 and 0.4999...
  p_temp = gt(abs(qfrac1_temp),1E-5); % indices of propagating phonons

  U_temp = LeadPhonon.U_minus;
  [qdata_temp, Udata_temp, vdata_temp, bdata_temp] = ...
    GenerateImageBZPoints(w,vlong_temp,qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,U_temp,Lead,LeadPhonDisp); % DEBUG

  LeadPhonon.U_minus = Udata_temp;
  LeadPhonon.VecQx_minus = qdata_temp(:,1);
  LeadPhonon.VecQy_minus = qdata_temp(:,2);
  LeadPhonon.VecQz_minus = qdata_temp(:,3);
  LeadPhonon.VecVelx_minus = vdata_temp(:,1);
  LeadPhonon.VecVely_minus = vdata_temp(:,2);
  LeadPhonon.VecVelz_minus = vdata_temp(:,3);
  LeadPhonon.VecB_minus  = bdata_temp(:,1);

  % =================================================================================================
  % Advanced minus phonon channels
  % =================================================================================================
  % === transverse wave vector(s) ===
  qfrac2_temp = [LeadPhonon.VecQ_tran1]/(2*pi); % fractional wave vector (transv. 1) (b/w 0 and 1)
  qfrac2_temp = mod(qfrac2_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...
  qfrac3_temp = [LeadPhonon.VecQ_tran2]/(2*pi); % fractional wave vector (transv. 2) (b/w 0 and 1)
  qfrac3_temp = mod(qfrac3_temp+0.5,1.0)-0.5;       % set value between -0.5 and 0.4999...

  % === longitudinal wave vector ===
  vlong_temp  = [LeadPhonon.VecV_minus_adv]; % AGF phonon longitudinal group velocity 
  qfrac1_temp = [LeadPhonon.VecQ_minus_adv]/(2*pi); % fractional wave vector (long.) (b/w 0 and 1)  
  qfrac1_temp = mod(qfrac1_temp+0.5,1.0)-0.5;           % set value between -0.5 and 0.4999...
  p_temp = gt(abs(qfrac1_temp),1E-5); % indices of propagating phonons

  U_temp = LeadPhonon.U_minus_adv;
  [qdata_temp, Udata_temp, vdata_temp, bdata_temp] = ...
    GenerateImageBZPoints(w,vlong_temp,qfrac1_temp,qfrac2_temp,qfrac3_temp,p_temp,U_temp,Lead,LeadPhonDisp); % DEBUG

  LeadPhonon.U_minus_adv = Udata_temp;
  LeadPhonon.VecQx_minus_adv = qdata_temp(:,1);
  LeadPhonon.VecQy_minus_adv = qdata_temp(:,2);
  LeadPhonon.VecQz_minus_adv = qdata_temp(:,3);
  LeadPhonon.VecVelx_minus_adv = vdata_temp(:,1);
  LeadPhonon.VecVely_minus_adv = vdata_temp(:,2);
  LeadPhonon.VecVelz_minus_adv = vdata_temp(:,3);
  LeadPhonon.VecB_minus_adv  = bdata_temp(:,1);

  % ------------------------------------------------------------------------------

  % LeadPhonon = rmfield(LeadPhonon,'VecQ_plus');
  % LeadPhonon = rmfield(LeadPhonon,'VecQ_plus_adv');
  % LeadPhonon = rmfield(LeadPhonon,'VecQ_minus');
  % LeadPhonon = rmfield(LeadPhonon,'VecQ_minus_adv');

  % =================================================================================================
  % Recompute the pristine transmission coeffs 
  % =================================================================================================
  t_L = 2*1i*w*sqrt(LeadPhonon.V_plus)*pinv(LeadPhonon.U_plus)*LeadPhonon.MatBulkG*pinv(LeadPhonon.U_minus_adv')*sqrt(LeadPhonon.V_minus_adv)/Lead.a_long;
  T_plus = diag(t_L*t_L');      % Perfect left-to-right transmission (out-going) in lead 
  T_minus_adv = diag(t_L'*t_L); % Perfect left-to-right transmission (in-coming) in lead 
  
  t_R = 2*1i*w*sqrt(LeadPhonon.V_minus)*pinv(LeadPhonon.U_minus)*LeadPhonon.MatBulkG*pinv(LeadPhonon.U_plus_adv')*sqrt(LeadPhonon.V_plus_adv)/Lead.a_long;
  T_minus = diag(t_R*t_R');     % Perfect right-to-left transmission (out-going) in lead 
  T_plus_adv = diag(t_R'*t_R);  % Perfect right-to-left transmission (in-coming) in lead 

  LeadPhonon.VecT_plus      = T_plus;
  LeadPhonon.VecT_plus_adv  = T_plus_adv;
  LeadPhonon.VecT_minus     = T_minus;
  LeadPhonon.VecT_minus_adv = T_minus_adv;
  LeadPhonon.Xi_mode = sum(T_plus);
end


