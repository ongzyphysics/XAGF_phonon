function ExtendedAtomisticGreensFunctionTransmission(InputFilesDir,OutputFilesDir)
% NOTE: InputFilesDir is the folder containing the input *.agf files
% NOTE: OutputFilesDir is the folder containing the output *.dat files 

  CurrDir = pwd; % current directory

  tic;

  % -----------------------------------------------------------------
  % Read input frequency data points
  % -----------------------------------------------------------------
  [wvec, wwvec] = ReadTransmissionFrequency(InputFilesDir);
  fprintf(1,'++ Input frequency range file <Input_Frequency.agf> read in %f seconds. \n\n',toc);
  % -----------------------------------------------------------------

  % -----------------------------------------------------------------
  % Read scattering results frequency data points
  % -----------------------------------------------------------------
  wvec_scatter = ReadScatteringFrequency(InputFilesDir);
  fprintf(1,'++ Scattering frequency range file <Scattering_Output_Frequency.agf> read in %f seconds. \n\n',toc);
  % -----------------------------------------------------------------

  % -----------------------------------------------------------------
  % Read force constant matrices, mass matrices and lattice constants 
  % -----------------------------------------------------------------
  tic;
  fprintf(1,'++ Input matrix files: \n');
  [LeftLeadParam, RightLeadParam, CenterParam] = ReadSimulationParameters(InputFilesDir); 
  fprintf(1,'   read in %f seconds. \n\n',toc);
  % -----------------------------------------------------------------

  % -----------------------------------------------------------------
  % Read phonon dispersion parameters and matrices 
  % -----------------------------------------------------------------
  tic;
  fprintf(1,'++ Input phonon dispersion files: \n');
  [LeftPhonDisp, RightPhonDisp] = ReadPhononDispersionParameters(InputFilesDir); 
  [LeftPhonDisp, RightPhonDisp] = ReadPhononMappingParameters(InputFilesDir,LeftPhonDisp,RightPhonDisp);
  fprintf(1,'   read in %f seconds. \n\n',toc);

  % -----------------------------------------------------------------
  % Convert force constants and atomic masses to H matrices (mass normalized matrices)
  % -----------------------------------------------------------------
  tic;
  [Left, Center, Right] = ConvertMassNormalizedMatrices(LeftLeadParam,CenterParam,RightLeadParam);
  fprintf(1,'++ Conversion of parameter took %f seconds. \n\n',toc);
  % -----------------------------------------------------------------

  for nw = 1:numel(wvec)
      % disp(nw);
      fprintf(1,'COMPUTATION: %4d of %4d at frequency = %10.4e rad/s \n\n',nw,length(wvec),wvec(nw));
      % -----------------------------------------------------------------
      % Calculate surface Green's functions for left lead
      % -----------------------------------------------------------------
      tic;
      if and(ge(nw,LeftLeadParam.nwread_min),le(nw,LeftLeadParam.nwread_max))
          fname_read = sprintf('%s/LeftLead_StoredData_%d.mat',InputFilesDir,nw);
          LeftPhonon(nw) = ReadStoredGFData(fname_read,wvec(nw));
          clear fname_read;
          fprintf(1,'++ Left surface Greens function loading took %f seconds. \n\n',toc); 
      else
          % LeftPhonon(nw) = ComputeSurfaceGreensFunction(wvec(nw),wwvec(nw),Left);
          TempLeftPhonon = ComputeSurfaceGreensFunction(wvec(nw),wwvec(nw),Left);  
            % This requires "TransverseSubspaceTransform.m" and "CombineTransverseSubspaceComponents.m"         
          % LeftPhonon(nw) = MapBulkPhononModes(wvec(nw),LeftPhonon(nw),Left,LeftPhonDisp);
          LeftPhonon(nw) = MapBulkPhononModes(wvec(nw),TempLeftPhonon,Left,LeftPhonDisp);
            % This requires "GenerateImageBZPoints.m"
          % LeftPhonon(nw) = RebuildExtendedModeEigenvector(wvec(nw),LeftPhonon(nw),Left,LeftPhonDisp);
          % LeftPhonon(nw) = SetPropagatingModeReflectionSymmetry(LeftPhonon(nw));
          fprintf(1,'++ Left surface Greens function computation took %f seconds. \n\n',toc); 
      end
           
      % -----------------------------------------------------------------

      % -----------------------------------------------------------------
      % Calculate surface Green's functions for right lead
      % -----------------------------------------------------------------      
      tic;
      if and(ge(nw,RightLeadParam.nwread_min),le(nw,RightLeadParam.nwread_max))
          fname_read = sprintf('%s/RightLead_StoredData_%d.mat',InputFilesDir,nw);
          RightPhonon(nw) = ReadStoredGFData(fname_read,wvec(nw));
          clear fname_read;
          fprintf(1,'++ Right surface Greens function loading took %f seconds. \n\n',toc);
      else          
          % RightPhonon(nw) = ComputeSurfaceGreensFunction(wvec(nw),wwvec(nw),Right);
          TempRightPhonon = ComputeSurfaceGreensFunction(wvec(nw),wwvec(nw),Right);
          % This requires "TransverseSubspaceTransform.m" and "CombineTransverseSubspaceComponents.m"
          % RightPhonon(nw) = MapBulkPhononModes(wvec(nw),RightPhonon(nw),Right,RightPhonDisp);
          RightPhonon(nw) = MapBulkPhononModes(wvec(nw),TempRightPhonon,Right,RightPhonDisp);
          % This requires "GenerateImageBZPoints.m"
          % RightPhonon(nw) = RebuildExtendedModeEigenvector(wvec(nw),RightPhonon(nw),Right,RightPhonDisp);
          % RightPhonon(nw) = SetPropagatingModeReflectionSymmetry(RightPhonon(nw));
          fprintf(1,'++ Right surface Greens function computation took %f seconds. \n\n',toc);
      end
      
      % -----------------------------------------------------------------

      % -----------------------------------------------------------------
      % Compute frequency-dependent transmission and reflection matrices
      % -----------------------------------------------------------------
      tic;
      PhononData(nw) = ComputeCenterPhononTransmission(wvec(nw),wwvec(nw),...
                       LeftPhonon(nw),RightPhonon(nw),Center,Left,Right); % This requires "GreensFunctionSolver.m"
      fprintf(1,'++ Phonon transmission computation took %f seconds. \n\n',toc);
      % -----------------------------------------------------------------

      % -----------------------------------------------------------------
      % Compute frequency-dependent phonon density of states for leads
      % -----------------------------------------------------------------
      tic;
      LeftPhononRho(nw) = GetPhononDensityOfStates(wvec(nw),LeftPhonon(nw));
      RightPhononRho(nw) = GetPhononDensityOfStates(wvec(nw),RightPhonon(nw));
      fprintf(1,'++ Phonon density of states computation took %f seconds. \n\n',toc);
      % -----------------------------------------------------------------

      % -----------------------------------------------------------------
      % Store GF data if required 
      % -----------------------------------------------------------------      
      if and(ge(nw,LeftLeadParam.nwdump_min),le(nw,LeftLeadParam.nwdump_max))
          tic;
          fname_dump = sprintf('%s/LeftLead_StoredData_%d.mat',InputFilesDir,nw);
          WriteStoredGFData(fname_dump,LeftPhonon(nw),wvec(nw));
          fprintf(1,'++ Left phonon GF data storage took %f seconds. \n\n',toc);
          clear Phonon_stored w_stored fname_dump;
      end

      if and(ge(nw,RightLeadParam.nwdump_min),le(nw,RightLeadParam.nwdump_max))
          tic;
          fname_dump = sprintf('%s/RightLead_StoredData_%d.mat',InputFilesDir,nw);
          WriteStoredGFData(fname_dump,RightPhonon(nw),wvec(nw));
          fprintf(1,'++ Right phonon GF data storage took %f seconds. \n\n',toc);
          clear Phonon_stored w_stored fname_dump;          
      end
      
      % -----------------------------------------------------------------
      % Remove data used for frequency-dependent transmission to save memory 
      % -----------------------------------------------------------------
      LeftPhonon(nw).MatSurfGL = [];
      LeftPhonon(nw).MatSurfGR = [];
      LeftPhonon(nw).MatBulkG  = [];
      LeftPhonon(nw).U_plus      = [];
      LeftPhonon(nw).U_plus_adv  = [];
      LeftPhonon(nw).U_minus     = [];
      LeftPhonon(nw).U_minus_adv = [];
      LeftPhonon(nw).V_plus      = [];
      LeftPhonon(nw).V_plus_adv  = [];
      LeftPhonon(nw).V_minus     = [];
      LeftPhonon(nw).V_minus_adv = [];

      RightPhonon(nw).MatSurfGL = [];
      RightPhonon(nw).MatSurfGR = [];
      RightPhonon(nw).MatBulkG  = [];
      RightPhonon(nw).U_plus      = [];
      RightPhonon(nw).U_plus_adv  = [];
      RightPhonon(nw).U_minus     = [];
      RightPhonon(nw).U_minus_adv = [];
      RightPhonon(nw).V_plus      = [];
      RightPhonon(nw).V_plus_adv  = [];
      RightPhonon(nw).V_minus     = [];
      RightPhonon(nw).V_minus_adv = [];
  end

  % save('TempData.mat','PhononData','LeftPhonon','RightPhonon','Left','Right','LeftPhonDisp','RightPhonDisp','wvec'); % DEBUG
  % return; % DEBUG LINE

  LeftPhonon = rmfield(LeftPhonon,{'MatSurfGL','MatSurfGR','MatBulkG', ...      
               'U_plus','U_plus_adv','U_minus','U_minus_adv', ... 
               'V_plus','V_plus_adv','V_minus','V_minus_adv'});
  RightPhonon = rmfield(RightPhonon,{'MatSurfGL','MatSurfGR','MatBulkG', ...      
               'U_plus','U_plus_adv','U_minus','U_minus_adv', ... 
               'V_plus','V_plus_adv','V_minus','V_minus_adv'});

  % -----------------------------------------------------------------
  % Postprocessing transmission data for left and right lead mode
  % -----------------------------------------------------------------
  fprintf(1,'POSTPROCESSING DATA: \n\n',nw,length(wvec));

  % -----------------------------------------------------------------
  % Write text files containing overall and individual phonon mode transmission results
  % Write text files containing phonon dispersion of left and right leads
  % -----------------------------------------------------------------
  tic;
  fprintf(1,'++ Output files: \n');
  WritePhononTransmissionOutput(OutputFilesDir,wvec,LeftPhonon,RightPhonon,PhononData);
  % WritePhononDispersionOutput(OutputFilesDir,Left,Right);
  WritePhononScatteringOutput(OutputFilesDir,wvec_scatter,wvec,LeftPhonon,RightPhonon,PhononData);
  WritePhononDensityOfStatesOutput(OutputFilesDir,wvec,LeftPhononRho,RightPhononRho);
  fprintf(1,'   saved in %f seconds. \n\n',toc);
  % -----------------------------------------------------------------
  disp(' ');

  return; % DEBUG LINE
end

