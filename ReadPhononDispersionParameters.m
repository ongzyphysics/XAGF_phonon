function [LeftPhonDisp, RightPhonDisp] = ReadPhononDispersionParameters(DataFilesDir)
% NOTE: Read from Left_Phonon_Parameters.agf and Right_Phonon_Parameters.agf 
% NOTE: Read primitive lattice vectors
% NOTE: Read force constant matrix information
% NOTE: Read Left_Phonon_K${n}.agf for n = 1,2,... and Left_Phonon_M.agf
% NOTE: Read Right_Phonon_K${m}.agf for m = 1,2,... and Right_Phonon_M.agf

  CurrDir = pwd; % current directory where the AGF code .m files are 

  Angstrom = 1E-10; % definition of angstrom in SI units (meter)

  cd(DataFilesDir); % change to AGF input directory

  % === Read left lead phonon parameters ===
  clear Rvec1 Rvec2 Rvec3 Gvec1 Gvec2 Gvec3 PrimCell 

  filename = 'Left_Phonon_Parameters.agf';
  fid = fopen(filename,'r');
  
  if gt(fid,0)
      fprintf(1,'          <%s> \n',filename);

      ScanData = textscan(fid,'%s','delimiter','\n','commentstyle','#');
      TextLineData = ScanData{1}; % cells of text line data with each cell corresponding to one text line
      clear ScanData;
      fclose(fid);
      
      % === Find positions of data BEGIN and END headers ===
      npos_begin_bz_lattice = 0; % line number for BEGIN_BZ_LATTICE
      npos_end_bz_lattice = 0; % line number for END_BZ_LATTICE      
      npos_begin_fc_matrix_index = 0; % line number for BEGIN_FC_MATRIX_INDEX
      npos_end_fc_matrix_index = 0; % line number for END_FC_MATRIX_INDEX       
      
      for nline = 1:numel(TextLineData)
          InputLine = TextLineData{nline}; % string containing characters of text line
          InputLine = [InputLine repmat(' ',1,50)]; % padded input text line
          
          CheckPhrase = 'BEGIN_BZ_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_bz_lattice = nline;
          end
          
          CheckPhrase = 'END_BZ_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_bz_lattice = nline;
          end
          
          CheckPhrase = 'BEGIN_FC_MATRIX_INDEX';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_fc_matrix_index = nline;
          end          
          
          CheckPhrase = 'END_FC_MATRIX_INDEX';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_fc_matrix_index = nline;
          end          
      end      
      
      % === Read primitive lattice vectors ===
      if and(npos_begin_bz_lattice>0,(npos_end_bz_lattice-npos_begin_bz_lattice)==4)
          InputLine = TextLineData{npos_begin_bz_lattice+1}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          Rvec1 = double([InputText{1} InputText{2} InputText{3}])*Angstrom; % primitive lattice vector in direction 1
          
          InputLine = TextLineData{npos_begin_bz_lattice+2}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          Rvec2 = double([InputText{1} InputText{2} InputText{3}])*Angstrom; % primitive lattice vector in direction 2

          InputLine = TextLineData{npos_begin_bz_lattice+3}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          Rvec3 = double([InputText{1} InputText{2} InputText{3}])*Angstrom; % primitive lattice vector in direction 3
      end
      
      % === Read force constant matrix information ===
      if and(npos_begin_fc_matrix_index>0,(npos_end_fc_matrix_index-npos_begin_fc_matrix_index)>=2)
          for nline = (npos_begin_fc_matrix_index+1):(npos_end_fc_matrix_index-1)
              InputLine = TextLineData{nline}; % read text line              
              InputText = textscan(InputLine,'%d %f %f %f ','CommentStyle','#');
              ncell = InputText{1};  % primitive unit cell index
              rfrac1 = InputText{2}; % fractional coordinate in lattice vector direction 1
              rfrac2 = InputText{3}; % fractional coordinate in lattice vector direction 2
              rfrac3 = InputText{4}; % fractional coordinate in lattice vector direction 3
              PrimCell(ncell).rvec = rfrac1*Rvec1 ...
                + rfrac2*Rvec2 + rfrac3*Rvec3; % lattice position of primitive unit cell
          end          
      end

      % === Read mass-normalized force constant matrices ===
      M = importdata('Left_Phonon_M.agf',' ',0);
      InvSqrtM = diag(1./sqrt(diag(M)));
      for ncell = 1:1:numel(PrimCell)
          fname = sprintf('Left_Phonon_K%d.agf',ncell); 
          TempK = importdata(fname,' ',0);
          nmatrows = size(TempK,1); % no. of rows
          nmatcols = size(TempK,2); % no. of columns

          if eq(nmatcols,2*nmatrows)
              K = TempK(1:nmatrows,1:nmatrows) ...
                  + 1i*TempK(1:nmatrows,(1:nmatrows)+nmatrows);
          else
              K = TempK;
          end
          PrimCell(ncell).H = InvSqrtM*K*InvSqrtM; % mass-normalized force-constant matrix for unit cell
      end

      Gvec1 = cross(Rvec2,Rvec3)/(Rvec1*cross(Rvec2,Rvec3)'); % primitive BZ reciprocal lattice vector 1
      Gvec2 = cross(Rvec3,Rvec1)/(Rvec2*cross(Rvec3,Rvec1)'); % primitive BZ reciprocal lattice vector 2
      Gvec3 = cross(Rvec1,Rvec2)/(Rvec3*cross(Rvec1,Rvec2)'); % primitive BZ reciprocal lattice vector 3

      LeftPhonDisp.PrimCell = PrimCell;
      LeftPhonDisp.Rvec1 = Rvec1;
      LeftPhonDisp.Rvec2 = Rvec2;
      LeftPhonDisp.Rvec3 = Rvec3;
      LeftPhonDisp.Gvec1 = Gvec1;
      LeftPhonDisp.Gvec2 = Gvec2;
      LeftPhonDisp.Gvec3 = Gvec3;      
  else
      error(sprintf('<!!!> %s not found.',filename));      
  end
  
  % === Read right lead phonon parameters ===
  clear Rvec1 Rvec2 Rvec3 Gvec1 Gvec2 Gvec3 PrimCell 

  filename = 'Right_Phonon_Parameters.agf';
  fid = fopen(filename,'r');
  
  if gt(fid,0)
      fprintf(1,'          <%s> \n',filename);

      ScanData = textscan(fid,'%s','delimiter','\n','commentstyle','#');
      TextLineData = ScanData{1}; % cells of text line data with each cell corresponding to one text line
      clear ScanData;
      fclose(fid);
      
      % === Find positions of data BEGIN and END headers ===
      npos_begin_bz_lattice = 0; % line number for BEGIN_BZ_LATTICE
      npos_end_bz_lattice = 0; % line number for END_BZ_LATTICE      
      npos_begin_fc_matrix_index = 0; % line number for BEGIN_FC_MATRIX_INDEX
      npos_end_fc_matrix_index = 0; % line number for END_FC_MATRIX_INDEX      
      
      for nline = 1:numel(TextLineData)
          InputLine = TextLineData{nline};
          InputLine = [InputLine repmat(' ',1,50)]; % padded input text line
          
          CheckPhrase = 'BEGIN_BZ_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_bz_lattice = nline;
          end
          
          CheckPhrase = 'END_BZ_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_bz_lattice = nline;
          end
          
          CheckPhrase = 'BEGIN_FC_MATRIX_INDEX';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_fc_matrix_index = nline;
          end          
          
          CheckPhrase = 'END_FC_MATRIX_INDEX';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_fc_matrix_index = nline;
          end          
      end      
      
      % === Read primitive lattice vectors ===
      if and(npos_begin_bz_lattice>0,(npos_end_bz_lattice-npos_begin_bz_lattice)==4)
          InputLine = TextLineData{npos_begin_bz_lattice+1}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          Rvec1 = double([InputText{1} InputText{2} InputText{3}])*Angstrom; % primitive lattice vector in direction 1
          
          InputLine = TextLineData{npos_begin_bz_lattice+2}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          Rvec2 = double([InputText{1} InputText{2} InputText{3}])*Angstrom; % primitive lattice vector in direction 2

          InputLine = TextLineData{npos_begin_bz_lattice+3}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          Rvec3 = double([InputText{1} InputText{2} InputText{3}])*Angstrom; % primitive lattice vector in direction 3
      end
      
      % === Read force constant matrix information ===
      if and(npos_begin_fc_matrix_index>0,(npos_end_fc_matrix_index-npos_begin_fc_matrix_index)>=2)
          for nline = (npos_begin_fc_matrix_index+1):(npos_end_fc_matrix_index-1)
              InputLine = TextLineData{nline}; % read text line              
              InputText = textscan(InputLine,'%d %f %f %f ','CommentStyle','#');
              ncell = InputText{1};  % primitive unit cell index
              rfrac1 = InputText{2}; % fractional coordinate in lattice vector direction 1
              rfrac2 = InputText{3}; % fractional coordinate in lattice vector direction 2
              rfrac3 = InputText{4}; % fractional coordinate in lattice vector direction 3
              PrimCell(ncell).rvec = rfrac1*Rvec1 ...
                + rfrac2*Rvec2 + rfrac3*Rvec3; % lattice position of primitive unit cell
          end          
      end

      % === Read mass-normalized force constant matrices ===
      M = importdata('Right_Phonon_M.agf',' ',0);
      InvSqrtM = diag(1./sqrt(diag(M)));
      for ncell = 1:1:numel(PrimCell)
          fname = sprintf('Right_Phonon_K%d.agf',ncell); 
          TempK = importdata(fname,' ',0);
          nmatrows = size(TempK,1); % no. of rows
          nmatcols = size(TempK,2); % no. of columns

          if eq(nmatcols,2*nmatrows)
              K = TempK(1:nmatrows,1:nmatrows) ...
                  + 1i*TempK(1:nmatrows,(1:nmatrows)+nmatrows);
          else
              K = TempK;
          end
          PrimCell(ncell).H = InvSqrtM*K*InvSqrtM; % mass-normalized force-constant matrix for unit cell
      end

      Gvec1 = cross(Rvec2,Rvec3)/(Rvec1*cross(Rvec2,Rvec3)'); % primitive BZ reciprocal lattice vector 1
      Gvec2 = cross(Rvec3,Rvec1)/(Rvec2*cross(Rvec3,Rvec1)'); % primitive BZ reciprocal lattice vector 2
      Gvec3 = cross(Rvec1,Rvec2)/(Rvec3*cross(Rvec1,Rvec2)'); % primitive BZ reciprocal lattice vector 3

      RightPhonDisp.PrimCell = PrimCell;
      RightPhonDisp.Rvec1 = Rvec1;
      RightPhonDisp.Rvec2 = Rvec2;
      RightPhonDisp.Rvec3 = Rvec3;
      RightPhonDisp.Gvec1 = Gvec1;
      RightPhonDisp.Gvec2 = Gvec2;
      RightPhonDisp.Gvec3 = Gvec3;      
  else
      error(sprintf('<!!!> %s not found.',filename));      
  end  

  cd(CurrDir); % return to directory containing AGF code .m files

end
