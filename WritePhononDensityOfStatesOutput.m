function WritePhononDensityOfStatesOutput(DataFilesDir,wvec,LeftPhononRho,RightPhononRho)
% NOTE: We write un-normalized left/right-lead phonon density of states data to "Output_Left_PhononDOS.dat" and "Output_Right_PhononDOS.dat"
% NOTE: Format of files: Frequency, density

  CurrDir = pwd; % current directory containing AGF code .m files

  filename_rho_L = 'Output_Left_PhononDOS.dat'; 
  filename_rho_R = 'Output_Right_PhononDOS.dat';

  cd(DataFilesDir); % change to AGF output directory
  fid_L = fopen(filename_rho_L,'w');
  fid_R = fopen(filename_rho_R,'w');
  cd(CurrDir); % return to directory containing AGF code .m files

  % disp(size(wvec));
  % disp(size(LeftPhononRho));
  % disp(size(RightPhononRho));

  nwmax = numel(wvec);

  if gt(nwmax,1)
      LeftNormFactor  = trapz(wvec,LeftPhononRho');
      RightNormFactor = trapz(wvec,RightPhononRho');

      for nw = 1:1:nwmax
          w = wvec(nw);
          fprintf(fid_L,'%12.4e %12.4e \n', w, LeftPhononRho(nw)/LeftNormFactor);
          fprintf(fid_R,'%12.4e %12.4e \n', w, RightPhononRho(nw)/RightNormFactor);
      end
  else
      for nw = 1:1:nwmax
          w = wvec(nw);
          fprintf(fid_L,'%12.4e %12.4e \n', w, 1.0/w);
          fprintf(fid_R,'%12.4e %12.4e \n', w, 1.0/w);
      end
  end

  fclose(fid_L);
  fclose(fid_R);

  fprintf(1,'\t  <%s> \n', filename_rho_L);
  fprintf(1,'\t  <%s> \n', filename_rho_R);
end
