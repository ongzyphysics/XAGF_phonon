function [LeftPhonDisp, RightPhonDisp] = ReadPhononMappingParameters(DataFilesDir,LeftPhonDisp_in,RightPhonDisp_in)
% NOTE: Read from Left_Phonon_Parameters.agf and Right_Phonon_Parameters.agf 
% NOTE: Read parameters for mapping of AGF lattice to primitive bulk lattice

  Angstrom = 1E-10;

  CurrDir = pwd; % current directory where the AGF code .m files are 
  cd(DataFilesDir); % change to AGF input directory

  LeftPhonDisp = LeftPhonDisp_in;
  RightPhonDisp = RightPhonDisp_in;

  LeftPhonDisp.nbasis_bulk = [];
  LeftPhonDisp.BasisAtomRxvec = [];
  LeftPhonDisp.BasisAtomRyvec = [];
  LeftPhonDisp.BasisAtomRzvec = [];

  RightPhonDisp.nbasis_bulk = [];
  RightPhonDisp.BasisAtomRxvec = [];
  RightPhonDisp.BasisAtomRyvec = [];
  RightPhonDisp.BasisAtomRzvec = [];

  % ---------------------------------------------------------------------  
  % Left bulk phonons
  % ---------------------------------------------------------------------
  filename = 'Left_Phonon_Parameters.agf';

  fid = fopen(filename,'r');

  if gt(fid,0)
      fprintf(1,'          <%s> \n',filename);

      ScanData = textscan(fid,'%s','delimiter','\n','commentstyle','#');
      TextLineData = ScanData{1}; % cells of text line data with each cell corresponding to one text line
      clear ScanData;
      fclose(fid);
      
      % === Find positions of data BEGIN and END headers ===
      npos_begin_agf2bulk_lattice = 0; % line number for BEGIN_AGF2BULK_LATTICE
      npos_end_agf2bulk_lattice = 0; % line number for END_AGF2BULK_LATTICE     
      
      for nline = 1:numel(TextLineData)
          InputLine = TextLineData{nline}; % string containing characters of text line
          InputLine = [InputLine repmat(' ',1,50)]; % padded input text line
          
          CheckPhrase = 'BEGIN_AGF2BULK_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_agf2bulk_lattice = nline;
          end
          
          CheckPhrase = 'END_AGF2BULK_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_agf2bulk_lattice = nline;
          end
      end      
      
      % === Read force constant matrix information ===
      if and(npos_begin_agf2bulk_lattice>0,(npos_end_agf2bulk_lattice-npos_begin_agf2bulk_lattice)>=2)
          for nline = (npos_begin_agf2bulk_lattice+1):(npos_end_agf2bulk_lattice-1)
              InputLine = TextLineData{nline}; % read text line                            
              InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
              nbasis_agf = InputText{1};  % AGF unit-cell basis atom index
              nbasis_bulk = InputText{2}; % Bulk lattice primitive unit-cell basis atom index
              rfrac1 = InputText{3}; % fractional coord of AGF unit cell wrt primitive lattice vector 1
              rfrac2 = InputText{4}; % fractional coord of AGF unit cell wrt primitive lattice vector 2
              rfrac3 = InputText{5}; % fractional coord of AGF unit cell wrt primitive lattice vector 3
              
              Rvec = rfrac1*LeftPhonDisp.Rvec1 + rfrac2*LeftPhonDisp.Rvec2 + rfrac3*LeftPhonDisp.Rvec3;
              LeftPhonDisp.nbasis_bulk(nbasis_agf) = nbasis_bulk;
              LeftPhonDisp.BasisAtomRxvec(nbasis_agf) = Rvec(1);
              LeftPhonDisp.BasisAtomRyvec(nbasis_agf) = Rvec(2);
              LeftPhonDisp.BasisAtomRzvec(nbasis_agf) = Rvec(3);
          end          
      end
  else
      error(sprintf('<!!!> %s not found.',filename));      
  end  

  % ---------------------------------------------------------------------  
  % Right bulk phonons
  % ---------------------------------------------------------------------
  filename = 'Right_Phonon_Parameters.agf';

  fid = fopen(filename,'r');
  
  if gt(fid,0)
      fprintf(1,'          <%s> \n',filename);

      ScanData = textscan(fid,'%s','delimiter','\n','commentstyle','#');
      TextLineData = ScanData{1}; % cells of text line data with each cell corresponding to one text line
      clear ScanData;
      fclose(fid);
      
      % === Find positions of data BEGIN and END headers ===
      npos_begin_agf2bulk_lattice = 0; % line number for BEGIN_AGF2BULK_LATTICE
      npos_end_agf2bulk_lattice = 0; % line number for END_AGF2BULK_LATTICE     
      
      for nline = 1:numel(TextLineData)
          InputLine = TextLineData{nline}; % string containing characters of text line
          InputLine = [InputLine repmat(' ',1,50)]; % padded input text line
          
          CheckPhrase = 'BEGIN_AGF2BULK_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_agf2bulk_lattice = nline;
          end
          
          CheckPhrase = 'END_AGF2BULK_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_agf2bulk_lattice = nline;
          end
      end      
      
      % === Read force constant matrix information ===
      if and(npos_begin_agf2bulk_lattice>0,(npos_end_agf2bulk_lattice-npos_begin_agf2bulk_lattice)>=2)
          for nline = (npos_begin_agf2bulk_lattice+1):(npos_end_agf2bulk_lattice-1)
              InputLine = TextLineData{nline}; % read text line                            
              InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
              nbasis_agf = InputText{1};  % AGF unit-cell basis atom index
              nbasis_bulk = InputText{2}; % Bulk lattice primitive unit-cell basis atom index
              rfrac1 = InputText{3}; % fractional coord of AGF unit cell wrt primitive lattice vector 1
              rfrac2 = InputText{4}; % fractional coord of AGF unit cell wrt primitive lattice vector 2
              rfrac3 = InputText{5}; % fractional coord of AGF unit cell wrt primitive lattice vector 3
              
              Rvec = rfrac1*RightPhonDisp.Rvec1 + rfrac2*RightPhonDisp.Rvec2 + rfrac3*RightPhonDisp.Rvec3;
              RightPhonDisp.nbasis_bulk(nbasis_agf) = nbasis_bulk;
              RightPhonDisp.BasisAtomRxvec(nbasis_agf) = Rvec(1);
              RightPhonDisp.BasisAtomRyvec(nbasis_agf) = Rvec(2);
              RightPhonDisp.BasisAtomRzvec(nbasis_agf) = Rvec(3);
          end          
      end
  else
      error(sprintf('<!!!> %s not found.',filename));      
  end  
  
  cd(CurrDir);
end
