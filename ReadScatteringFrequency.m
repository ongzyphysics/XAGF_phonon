function wvec_scatter = ReadScatteringFrequency(DataFilesDir)
% NOTE: Reads the frequency points for outputting channel scattering data
% NOTE: All scattering frequency points should coincide with some or all transmission frequency points
% NOTE: This is optional. If file is not found, then we assume no channel scattering data output is wanted
% NOTE: This tells us the frequency points where we calculate the S matrix

  CurrDir = pwd; % current directory containing AGF code .m files

  cd(DataFilesDir); % change to AGF input directory

  % -----------------------------------------------------------------
  % Read data file
  % -----------------------------------------------------------------
  filename = 'Scattering_Output_Frequency.agf'; % data file name containing scattering frequency data points
  fid = fopen(filename,'r');

  if gt(fid,0)
      data = importdata(filename,' ',0); % read file containing range of frequency values  
      wvec_scatter = data(:,1);
  else
      warning(sprintf('<!!!> %s not found.',filename));
      wvec_scatter = 0;
  end
  fclose(fid);

  cd(CurrDir); % return to AGF code directory
end
