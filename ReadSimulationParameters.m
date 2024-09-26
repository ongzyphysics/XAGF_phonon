function [LeftLeadParam, RightLeadParam, CenterParam] = ReadSimulationParameters(DataFilesDir)
% NOTE: Read KL, KC, KR and M matrices for left lead in Left_KL.agf, Left_KC.agf, Left_KR.agf and Left_M.agf 
% NOTE: Read KL, KC, KR and M matrices for right lead in Right_KL.agf, Right_KC.agf, Right_KR.agf and Right_M.agf 
% NOTE: Read KL, KC, KR and M matrices for each layer in center in Center_KL${n}.agf, Center_KC${n}.agf, 
% NOTE:   Center_KR${n}.agf and Center_M${n}.agf where n is the layer number index
% NOTE: Read Left_Parameters.agf and  Right_Parameters.agf for: 
% NOTE:   (1) non-primitive lattice vectors of AGF set up,
% NOTE:   (2) cross-sectional span of cells in transverse directions, and 
% NOTE:   (3) arrangement of cells in transverse directions

  angstrom = 1.0E-10; % definition of angstrom in SI units (meter)

  CurrDir = pwd; % current directory containing AGF code .m files
  cd(DataFilesDir); % change to AGF input directory

  % -------------------------------------------------------------------
  % Read parameters for left lead (KL, KC, KR & M matrices and parameters)
  %   The no. of rows = no. of columns if there are no imaginary components
  %   The no. of rows = 2 x no. of columns if there are imaginary components
  % -------------------------------------------------------------------

  % === Left-lead KL matrix ===
  filename = 'Left_KL.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'\t  <%s> ',filename); % name of text file containing KL matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          LeftLeadParam.MatKL = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatKL = TempMat;
      end    
      LeftLeadParam.MatKL = sparse(LeftLeadParam.MatKL);
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing KL matrix
      tmp = load(filename_mat,'-mat');
      LeftLeadParam.MatKL = sparse(tmp.KL);
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));
  end

  % === Left-lead KC matrix ===
  filename = 'Left_KC.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'   <%s> ',filename); % name of text file containing KC matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          LeftLeadParam.MatKC = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatKC = TempMat;
      end    
      LeftLeadParam.MatKC = sparse(LeftLeadParam.MatKC);
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing KC matrix
      tmp = load(filename_mat,'-mat');
      LeftLeadParam.MatKC = sparse(tmp.KC);
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));
  end

  % === Left-lead KR matrix ===
  filename = 'Left_KR.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'   <%s> ',filename); % name of text file containing KR matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          LeftLeadParam.MatKR = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatKR = TempMat;
      end
      LeftLeadParam.MatKR = sparse(LeftLeadParam.MatKR);    
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing KR matrix  
      tmp = load(filename_mat,'-mat');
      LeftLeadParam.MatKR = sparse(tmp.KR);
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));
  end

  % === Left-lead M matrix ===
  filename = 'Left_M.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'   <%s> ',filename);

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          LeftLeadParam.MatM = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          LeftLeadParam.MatM = TempMat;
      end    
      LeftLeadParam.MatM = sparse(LeftLeadParam.MatM);
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing M matrix
      tmp = load(filename_mat,'-mat');
      LeftLeadParam.MatM = sparse(tmp.M);
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));
  end

  % === Left lead parameters ===
  filename = 'Left_Parameters.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s> \n',filename);
      
      ScanData = textscan(fid,'%s','delimiter','\n','commentstyle','#');
      TextLineData = ScanData{1}; % cells of text line data with each cell corresponding to one text line
      clear ScanData;
      fclose(fid);
      
      % === Find positions of data BEGIN and END headers ===
      npos_begin_lattice = 0; % line number for BEGIN_LATTICE
      npos_end_lattice = 0; % line number for END_LATTICE
      npos_begin_transverse_span = 0; % line number for BEGIN_TRANSVERSE_SPAN
      npos_end_transverse_span = 0; % line number for END_TRANSVERSE_SPAN
      npos_begin_cell_index = 0; % line number for BEGIN_CELL_INDEX
      npos_end_cell_index = 0; % line number for END_CELL_INDEX
      npos_begin_gf_read_data = 0; % line number for BEGIN_GF_READ_DATA
      npos_end_gf_read_data = 0; % line number for END_GF_READ_DATA
      npos_begin_gf_dump_data = 0; % line number for BEGIN_GF_DUMP_DATA
      npos_end_gf_dump_data = 0; % line number for END_GF_DUMP_DATA
      
      for nline = 1:numel(TextLineData)
          InputLine = TextLineData{nline};
          InputLine = [InputLine repmat(' ',1,50)]; % padded input text line
          
          CheckPhrase = 'BEGIN_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_lattice = nline;
          end
          
          CheckPhrase = 'END_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_lattice = nline;
          end
          
          CheckPhrase = 'BEGIN_TRANSVERSE_SPAN';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_transverse_span = nline;
          end          
          
          CheckPhrase = 'END_TRANSVERSE_SPAN';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_transverse_span = nline;
          end          
          
          CheckPhrase = 'BEGIN_CELL_INDEX';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_cell_index = nline;
          end
          
          CheckPhrase = 'END_CELL_INDEX';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_cell_index = nline;
          end
          
          CheckPhrase = 'BEGIN_GF_READ_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_gf_read_data = nline;
          end          
          
          CheckPhrase = 'END_GF_READ_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_gf_read_data = nline;
          end          
          
          CheckPhrase = 'BEGIN_GF_DUMP_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_gf_dump_data = nline;
          end          
          
          CheckPhrase = 'END_GF_DUMP_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_gf_dump_data = nline;
          end          
      end
      
      % === Set lattice (non-primitive) vectors ===
      if and(npos_begin_lattice>0,(npos_end_lattice-npos_begin_lattice)==4)
          InputLine = TextLineData{npos_begin_lattice+1}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          LeftLeadParam.rvec_long = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          % lattice vector in longitudinal direction
          LeftLeadParam.a_long = norm(LeftLeadParam.rvec_long);
          % lattice constant in longitudinal direction
          
          InputLine = TextLineData{npos_begin_lattice+2}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          LeftLeadParam.rvec_tran1 = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          % lattice vector in transverse direction 1
          
          InputLine = TextLineData{npos_begin_lattice+3}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          LeftLeadParam.rvec_tran2 = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          % lattice vector in transverse direction 2          
      end
      
      % === Set tranverse span ===
      if and(npos_begin_transverse_span>0,(npos_end_transverse_span-npos_begin_transverse_span)==2)
          InputLine = TextLineData{npos_begin_transverse_span+1}; % read text line
          InputText = textscan(InputLine,'%d %d %d ','CommentStyle','#');
          LeftLeadParam.n_tran1 = InputText{2}; % no. of unit cells in transverse direction 1
          LeftLeadParam.n_tran2 = InputText{3}; % no. of unit cells in transverse direction 2
      end     
      
      % === Set tranverse cell information ===
      if and(npos_begin_cell_index>0,(npos_end_cell_index-npos_begin_cell_index)>=2)
          for nline = (npos_begin_cell_index+1):(npos_end_cell_index-1)
              InputLine = TextLineData{nline}; % read text line              
              InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
              ncell = InputText{1};  % transverse cell index
              natoms = InputText{2}; % number of atoms in transverse cell
              LeftLeadParam.natomdim = size(LeftLeadParam.MatM,1)/(LeftLeadParam.n_tran1*LeftLeadParam.n_tran2*natoms); % no. of DOFs per atom
              LeftLeadParam.TransCell(ncell).nrange = (1:(LeftLeadParam.natomdim*natoms)) + (ncell-1)*LeftLeadParam.natomdim*natoms;
              frac_long = InputText{3};  % fractional coordinate in longitudinal direction
              frac_tran1 = InputText{4}; % fractional coordinate in transverse direction 1
              frac_tran2 = InputText{5}; % fractional coordinate in transverse direction 2
              LeftLeadParam.TransCell(ncell).rvec = frac_long*LeftLeadParam.rvec_long ...
                  + frac_tran1*LeftLeadParam.rvec_tran1 + frac_tran2*LeftLeadParam.rvec_tran2;
              % lattice position of transverse unit cell              
          end          
      end
      
      % === Set pre-processed lead data information (read) ===
      if and(npos_begin_gf_read_data>0,(npos_end_gf_read_data-npos_begin_gf_read_data)==2) 
          InputLine = TextLineData{npos_begin_gf_read_data+1}; % read text line
          InputText = textscan(InputLine,'%d %d ','CommentStyle','#');
          LeftLeadParam.nwread_min = InputText{1}; % initial index of freq-depend. surface Green's functions
          LeftLeadParam.nwread_max = InputText{2}; % final index of freq-depend. surface Green's functions          
      end
      
      % === Set post-processed lead data information (write) ===
      if and(npos_begin_gf_dump_data>0,(npos_end_gf_dump_data-npos_begin_gf_dump_data)==2) 
          InputLine = TextLineData{npos_begin_gf_dump_data+1}; % read text line
          InputText = textscan(InputLine,'%d %d ','CommentStyle','#');
          LeftLeadParam.nwdump_min = InputText{1}; % initial index of freq-depend. surface Green's functions
          LeftLeadParam.nwdump_max = InputText{2}; % final index of freq-depend. surface Green's functions          
      end      
      
      clear TextLineData InputLine InputText;
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % -------------------------------------------------------------------
  % Read parameters for right lead (KL, KC, KR & M matrices and parameters)
  %   The no. of rows = no. of columns if there are no imaginary components
  %   The no. of rows = 2 x no. of columns if there are imaginary components
  % -------------------------------------------------------------------

  % === Right KL matrix ===
  filename = 'Right_KL.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'\t  <%s> ',filename); % name of text file containing KL matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          RightLeadParam.MatKL = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          RightLeadParam.MatKL = TempMat;
      end  
      RightLeadParam.MatKL = sparse(RightLeadParam.MatKL);              
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing KL matrix  
      tmp = load(filename_mat,'-mat');
      RightLeadParam.MatKL = sparse(tmp.KL);              
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));      
  end

  % === Right KC matrix ===
  filename = 'Right_KC.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'   <%s>',filename); % name of text file containing KC matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          RightLeadParam.MatKC = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          RightLeadParam.MatKC = TempMat;
      end    
      RightLeadParam.MatKC = sparse(RightLeadParam.MatKC);              
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing KC matrix    
      tmp = load(filename_mat,'-mat');
      RightLeadParam.MatKC = sparse(tmp.KC);              
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));      
  end

  % === Right KR matrix ===
  filename = 'Right_KR.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'   <%s>',filename); % name of text file containing KR matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          RightLeadParam.MatKR = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          RightLeadParam.MatKR = TempMat;
      end    
      RightLeadParam.MatKR = sparse(RightLeadParam.MatKR);    
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing KR matrix      
      tmp = load(filename_mat,'-mat');
      RightLeadParam.MatKR = sparse(tmp.KR);              
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));      
  end

  % === Right M matrix ===
  filename = 'Right_M.agf';
  filename_mat = sprintf('%sm',filename); % file name of matlab data (.agfm) version 
  if exist(filename)==2
      fprintf(1,'   <%s>',filename); % name of text file containing M matrix

      TempMat = importdata(filename,' ',0);
      nmatrows = size(TempMat,1); % no. of rows
      nmatcols = size(TempMat,2); % no. of columns

      if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
          RightLeadParam.MatM = ...
              TempMat(1:nmatrows,1:nmatrows) ...
              + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
      else % i.e. if there are no imaginary components
          RightLeadParam.MatM = TempMat;
      end    
      RightLeadParam.MatM = sparse(RightLeadParam.MatM);                    
  elseif exist(filename_mat)==2
      fprintf(1,'\t  <%s> ',filename_mat); % name of text file containing M matrix      
      tmp = load(filename_mat,'-mat');
      RightLeadParam.MatM = sparse(tmp.M);                    
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));      
  end

  % === Right lead parameters ===
  filename = 'Right_Parameters.agf';
  fid = fopen(filename,'r');
  if gt(fid,0)
      fprintf(1,'   <%s> \n',filename);    
      
      ScanData = textscan(fid,'%s','delimiter','\n','commentstyle','#');
      TextLineData = ScanData{1}; % cells of text line data with each cell corresponding to one text line
      clear ScanData;
      fclose(fid);
      
      % === Find positions of data BEGIN and END headers ===
      npos_begin_lattice = 0; % line number for BEGIN_LATTICE
      npos_end_lattice = 0; % line number for END_LATTICE
      npos_begin_transverse_span = 0; % line number for BEGIN_TRANSVERSE_SPAN
      npos_end_transverse_span = 0; % line number for END_TRANSVERSE_SPAN
      npos_begin_cell_index = 0; % line number for BEGIN_CELL_INDEX
      npos_end_cell_index = 0; % line number for END_CELL_INDEX
      npos_begin_gf_read_data = 0; % line number for BEGIN_GF_READ_DATA
      npos_end_gf_read_data = 0; % line number for END_GF_READ_DATA
      npos_begin_gf_dump_data = 0; % line number for BEGIN_GF_DUMP_DATA
      npos_end_gf_dump_data = 0; % line number for END_GF_DUMP_DATA
      
      for nline = 1:numel(TextLineData)
          InputLine = TextLineData{nline};
          InputLine = [InputLine repmat(' ',1,50)]; % padded input text line
          
          CheckPhrase = 'BEGIN_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_lattice = nline;
          end
          
          CheckPhrase = 'END_LATTICE';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_lattice = nline;
          end
          
          CheckPhrase = 'BEGIN_TRANSVERSE_SPAN';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_transverse_span = nline;
          end          
          
          CheckPhrase = 'END_TRANSVERSE_SPAN';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_transverse_span = nline;
          end          
          
          CheckPhrase = 'BEGIN_CELL_INDEX';
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_cell_index = nline;
          end
          
          CheckPhrase = 'END_CELL_INDEX';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_cell_index = nline;
          end
          
          CheckPhrase = 'BEGIN_GF_READ_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_gf_read_data = nline;
          end          
          
          CheckPhrase = 'END_GF_READ_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_gf_read_data = nline;
          end          
          
          CheckPhrase = 'BEGIN_GF_DUMP_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_begin_gf_dump_data = nline;
          end          
          
          CheckPhrase = 'END_GF_DUMP_DATA';          
          if strcmp(InputLine(1:numel(CheckPhrase)),CheckPhrase)
              npos_end_gf_dump_data = nline;
          end          
      end
      
      % === Set lattice (non-primitive) vectors ===
      if and(npos_begin_lattice>0,(npos_end_lattice-npos_begin_lattice)==4)
          InputLine = TextLineData{npos_begin_lattice+1}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          RightLeadParam.rvec_long = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          % lattice vector in longitudinal direction
          RightLeadParam.a_long = norm(RightLeadParam.rvec_long);
          % lattice constant in longitudinal direction
          
          InputLine = TextLineData{npos_begin_lattice+2}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          RightLeadParam.rvec_tran1 = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          % lattice vector in transverse direction 1
          
          InputLine = TextLineData{npos_begin_lattice+3}; % read text line
          InputText = textscan(InputLine,'%f %f %f ','CommentStyle','#');
          RightLeadParam.rvec_tran2 = double([InputText{1} InputText{2} InputText{3}])*angstrom;
          % lattice vector in transverse direction 2          
      end
      
      % === Set tranverse span ===
      if and(npos_begin_transverse_span>0,(npos_end_transverse_span-npos_begin_transverse_span)==2)
          InputLine = TextLineData{npos_begin_transverse_span+1}; % read text line
          InputText = textscan(InputLine,'%d %d %d ','CommentStyle','#');
          RightLeadParam.n_tran1 = InputText{2}; % no. of unit cells in transverse direction 1
          RightLeadParam.n_tran2 = InputText{3}; % no. of unit cells in transverse direction 2
      end     
      
      % === Set tranverse cell information ===
      if and(npos_begin_cell_index>0,(npos_end_cell_index-npos_begin_cell_index)>=2)
          for nline = (npos_begin_cell_index+1):(npos_end_cell_index-1)
              InputLine = TextLineData{nline}; % read text line              
              InputText = textscan(InputLine,'%d %d %f %f %f ','CommentStyle','#');
              ncell = InputText{1};  % transverse cell index
              natoms = InputText{2}; % number of atoms in transverse cell
              RightLeadParam.natomdim = size(RightLeadParam.MatM,1)/(RightLeadParam.n_tran1*RightLeadParam.n_tran2*natoms); % no. of DOFs per atom
              RightLeadParam.TransCell(ncell).nrange = (1:(RightLeadParam.natomdim*natoms)) + (ncell-1)*RightLeadParam.natomdim*natoms;
              frac_long = InputText{3};  % fractional coordinate in longitudinal direction
              frac_tran1 = InputText{4}; % fractional coordinate in transverse direction 1
              frac_tran2 = InputText{5}; % fractional coordinate in transverse direction 2
              RightLeadParam.TransCell(ncell).rvec = frac_long*RightLeadParam.rvec_long ...
                  + frac_tran1*RightLeadParam.rvec_tran1 + frac_tran2*RightLeadParam.rvec_tran2;
              % lattice position of transverse unit cell              
          end          
      end
      
      % === Set pre-processed lead data information (read) ===
      if and(npos_begin_gf_read_data>0,(npos_end_gf_read_data-npos_begin_gf_read_data)==2) 
          InputLine = TextLineData{npos_begin_gf_read_data+1}; % read text line
          InputText = textscan(InputLine,'%d %d ','CommentStyle','#');
          RightLeadParam.nwread_min = InputText{1}; % initial index of freq-depend. surface Green's functions
          RightLeadParam.nwread_max = InputText{2}; % final index of freq-depend. surface Green's functions          
      end
      
      % === Set post-processed lead data information (write) ===
      if and(npos_begin_gf_dump_data>0,(npos_end_gf_dump_data-npos_begin_gf_dump_data)==2) 
          InputLine = TextLineData{npos_begin_gf_dump_data+1}; % read text line
          InputText = textscan(InputLine,'%d %d ','CommentStyle','#');
          RightLeadParam.nwdump_min = InputText{1}; % initial index of freq-depend. surface Green's functions
          RightLeadParam.nwdump_max = InputText{2}; % final index of freq-depend. surface Green's functions          
      end      
      
      clear TextLineData InputLine InputText;
      % fclose(fid);
  else
      error(sprintf('<!!!> %s not found.',filename));
  end

  % -------------------------------------------------------------------
  % Read parameters for channel/scattering region 
  % (KL, KC, KR and M submatrices)  
  %
  %     To facilitate the use of the recursive Green's function (RCF) 
  %     algorithm, we partition, if possible, the scattering region 
  %     into successive layers with their own KL, KC, KR and M 
  %     submatrices. This partitioning  converts the giant matrices 
  %     associated with the scattering region into block-tridiagonal 
  %     matrices such that the diagonal and off-diagonal submatrices 
  %     correspond to the aforementioned KL, KC, KR and M matrices. The 
  %     partitioning also reduces the costs of the matrix inversion 
  %     operations.
  %
  % -------------------------------------------------------------------

  nlayer = 1;
  filename = sprintf('Center_M%d.agf',nlayer);
  filename_mat = sprintf('Center_M%d.agfm',nlayer);
  % fid = fopen(filename);

  if exist(filename)==2   
      while( exist(filename)==2 )
          % === M submatrix for layer 'nlayer' ===
          M_filename = sprintf('Center_M%d.agf',nlayer);

          TempMat = importdata(M_filename,' ',0);
          nmatrows = size(TempMat,1); % no. of rows
          nmatcols = size(TempMat,2); % no. of columns


          if eq(nmatcols,2*nmatrows) % i.e. if there are imaginary components
              Lyr(nlayer).MatM = ...
                  TempMat(1:nmatrows,1:nmatrows) ...
                  + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
          else
              Lyr(nlayer).MatM = TempMat;
          end    
          Lyr(nlayer).MatM = sparse(Lyr(nlayer).MatM);

          % fclose(fid); % close text file with M matrix data        
          nlayer = nlayer + 1;
        	
          filename = sprintf('Center_M%d.agf',nlayer);
          % fid = fopen(filename); % try to open next text file with M matrix data
      end
  elseif exist(filename_mat)==2 % Read M matrix data for layer 'nlayer'
      while exist(filename_mat)==2
          tmp = load(filename_mat,'-mat','M');
          % Lyr(nlayer).MatM = full(tmp.M);
          Lyr(nlayer).MatM = sparse(tmp.M);
          
          nlayer = nlayer + 1;
          filename_mat = sprintf('Center_M%d.agfm',nlayer);
      end
  else
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));
  end

  % -------------------------------------------------------------------

  nlayer = 1;
  filename = sprintf('Center_KC%d.agf',nlayer);
  filename_mat = sprintf('Center_KC%d.agfm',nlayer);
  % fid = fopen(filename);

  if exist(filename)==2
      while( exist(filename)==2 )
          if eq(nlayer,1)
              nrows_L = size(LeftLeadParam.MatM,1); % no. of degrees of freedom in slice left of current layer
          else
              nrows_L = size(Lyr(nlayer-1).MatM,1); % no. of degrees of freedom in slice left of current layer
          end

          nrows_C = size(Lyr(nlayer).MatM,1); % no. of degrees of freedom in current layer

          if eq(nlayer,numel(Lyr))
              nrows_R = size(RightLeadParam.MatM,1); % no. of degrees of freedom in slice right of current layer
          else
              nrows_R = size(Lyr(nlayer+1).MatM,1); % no. of degrees of freedom in slice right of current layer
          end

          % === KL submatrix for layer 'nlayer' ===
          KL_filename = sprintf('Center_KL%d.agf',nlayer);
          if not(exist(KL_filename)==2)
              error(sprintf('<!!!> %s not found.',KL_filename));
          end

          TempMat = importdata(KL_filename,' ',0);
          nmatrows = size(TempMat,1); % no. of rows
          nmatcols = size(TempMat,2); % no. of columns

          if eq(nmatcols,2*nrows_L)
              Lyr(nlayer).MatKL = ...
                  TempMat(1:nrows_C,1:nrows_L) ...
                  + 1i*TempMat(1:nrows_C,(1:nrows_L)+nrows_L);
          else
              Lyr(nlayer).MatKL = TempMat;
          end    
          Lyr(nlayer).MatKL = sparse(Lyr(nlayer).MatKL);

          % === KC submatrix for layer 'nlayer' ===
          KC_filename = sprintf('Center_KC%d.agf',nlayer);
          if not(exist(KC_filename)==2)
              error(sprintf('<!!!> %s not found.',KC_filename));
          end

          TempMat = importdata(KC_filename,' ',0);
          nmatrows = size(TempMat,1); % no. of rows
          nmatcols = size(TempMat,2); % no. of columns

          if eq(nmatcols,2*nrows_C)
              Lyr(nlayer).MatKC = ...
                  TempMat(1:nrows_C,1:nrows_C) ...
                  + 1i*TempMat(1:nrows_C,(1:nrows_C)+nrows_C);
          else
              Lyr(nlayer).MatKC = TempMat;
          end  
          Lyr(nlayer).MatKC = sparse(Lyr(nlayer).MatKC);  

          % === KR submatrix for layer 'nlayer' ===
          KR_filename = sprintf('Center_KR%d.agf',nlayer);
          if not(exist(KR_filename)==2)
              error(sprintf('<!!!> %s not found.',KR_filename));
          end         

          TempMat = importdata(KR_filename,' ',0);
          nmatrows = size(TempMat,1); % no. of rows
          nmatcols = size(TempMat,2); % no. of columns

          if eq(nmatcols,2*nrows_R)
              Lyr(nlayer).MatKR = ...
                  TempMat(1:nrows_C,1:nrows_R) ...
                  + 1i*TempMat(1:nrows_C,(1:nrows_R)+nrows_R);
          else
              Lyr(nlayer).MatKR = TempMat;
          end    
          Lyr(nlayer).MatKR = sparse(Lyr(nlayer).MatKR);

          % === M submatrix for layer 'nlayer' ===
          M_filename = sprintf('Center_M%d.agf',nlayer);
          if not(exist(M_filename)==2)
              error(sprintf('<!!!> %s not found.',M_filename));
          end          

          TempMat = importdata(M_filename,' ',0);
          nmatrows = size(TempMat,1); % no. of rows
          nmatcols = size(TempMat,2); % no. of columns

          if eq(nmatcols,2*nmatrows)
              Lyr(nlayer).MatM = ...
                  TempMat(1:nmatrows,1:nmatrows) ...
                  + 1i*TempMat(1:nmatrows,(1:nmatrows)+nmatrows);
          else
              Lyr(nlayer).MatM = TempMat;
          end    
          Lyr(nlayer).MatM = sparse(Lyr(nlayer).MatM);

          fprintf(1,'\t  <%s>', KL_filename);
          fprintf(1,'   <%s>', KC_filename);
          fprintf(1,'   <%s>', KR_filename);
          fprintf(1,'   <%s> \n', M_filename);            
        
          % fclose(fid);        
          nlayer = nlayer + 1;
        
          filename = sprintf('Center_KC%d.agf',nlayer);
          % fid = fopen(filename);
      end
  elseif exist(filename_mat)==2      
      while( exist(filename_mat)==2 )  
          % === KL submatrix for layer 'nlayer' ===
          KL_filename_mat = sprintf('Center_KL%d.agfm',nlayer);
          if not(exist(KL_filename_mat)==2)
              error(sprintf('<!!!> %s not found.',KL_filename_mat));
          end                    

          tmp = load(KL_filename_mat,'-mat','KL');
          % Lyr(nlayer).MatKL = full(tmp.KL);
          Lyr(nlayer).MatKL = sparse(tmp.KL);

          % === KC submatrix for layer 'nlayer' ===
          KC_filename_mat = sprintf('Center_KC%d.agfm',nlayer);
          if not(exist(KC_filename_mat)==2)
              error(sprintf('<!!!> %s not found.',KC_filename_mat));
          end                    

          tmp = load(KC_filename_mat,'-mat','KC');
          % Lyr(nlayer).MatKC = full(tmp.KC);
          Lyr(nlayer).MatKC = sparse(tmp.KC);

          % === KR submatrix for layer 'nlayer' ===
          KR_filename_mat = sprintf('Center_KR%d.agfm',nlayer);
          if not(exist(KR_filename_mat)==2)
              error(sprintf('<!!!> %s not found.',KR_filename_mat));
          end                              

          tmp = load(KR_filename_mat,'-mat','KR');
          % Lyr(nlayer).MatKR = full(tmp.KR);
          Lyr(nlayer).MatKR = sparse(tmp.KR);

          % === M submatrix for layer 'nlayer' ===
          M_filename_mat = sprintf('Center_M%d.agfm',nlayer);
          if not(exist(M_filename_mat)==2)
              error(sprintf('<!!!> %s not found.',M_filename_mat));
          end                              

          tmp = load(M_filename_mat,'-mat','M');
          % Lyr(nlayer).MatM = full(tmp.M);
          Lyr(nlayer).MatM = sparse(tmp.M);

          fprintf(1,'\t  <%s>', KL_filename_mat);
          fprintf(1,'   <%s>', KC_filename_mat);
          fprintf(1,'   <%s>', KR_filename_mat);
          fprintf(1,'   <%s> \n', M_filename_mat);            
          
          nlayer = nlayer + 1; 
          filename_mat = sprintf('Center_KC%d.agfm',nlayer);
      end
  else 
      error(sprintf('<!!!> %s and %s not found.',filename,filename_mat));
  end

  for nlayer = 1:numel(Lyr) % DEBUG
      % disp(nlayer)
      % disp(size(Lyr(nlayer).MatKL));
      % disp(size(Lyr(nlayer).MatKC));
      % disp(size(Lyr(nlayer).MatKR));
  end

  % -------------------------------------------------------------------
  % Padding of additional layers to the left and right edges of the 
  % scattering region. 
  % -------------------------------------------------------------------

  nlayer= numel(Lyr); % original no. of principal layers, this will be increased by 2 because of padding

  for n = (nlayer+1):-1:2 % move layer 1 to 2, 2 to 3, ..., nlayer to nlayer+1
      Lyr(n) = Lyr(n-1);
  end 

  % Set layer 1 (left edge) as layer from left lead
  Lyr(1).MatM  = LeftLeadParam.MatM;
  Lyr(1).MatKL = LeftLeadParam.MatKL;
  Lyr(1).MatKC = LeftLeadParam.MatKC;
  Lyr(1).MatKR = (Lyr(2).MatKL)';

  % Set layer nlayer+2 (right edge) as layer from right lead
  Lyr(nlayer+2).MatM  = RightLeadParam.MatM;
  Lyr(nlayer+2).MatKL = (Lyr(nlayer+1).MatKR)';
  Lyr(nlayer+2).MatKC = RightLeadParam.MatKC;
  Lyr(nlayer+2).MatKR = RightLeadParam.MatKR;

  nlayer= numel(Lyr); 

  CenterParam.nlayers = length(Lyr);
  CenterParam.Lyr = Lyr;
  CenterParam.MatKCL = Lyr(1).MatKL;
  CenterParam.MatKLC = (CenterParam.MatKCL)';
  CenterParam.MatKCR = Lyr(nlayer).MatKR;
  CenterParam.MatKRC = (CenterParam.MatKCR)';

  % -------------------------------------------------------------------
  % Message on number of degrees per atom
  % -------------------------------------------------------------------
  fprintf(1,'\t  Determined %d degrees of freedom per atom in left lead. \n', LeftLeadParam.natomdim);
  fprintf(1,'\t  Determined %d degrees of freedom per atom in right lead. \n', RightLeadParam.natomdim);  

  % -------------------------------------------------------------------
  % Message on GF data read and dump
  % -------------------------------------------------------------------  
  if ~and(isfield(LeftLeadParam,'nwread_min'),isfield(LeftLeadParam,'nwread_max'))
      LeftLeadParam.nwread_min = 0;
      LeftLeadParam.nwread_max = 0;
  end
  fprintf(1,'\t  To read left-lead GF data for frequency points: %4d to %4d \n', LeftLeadParam.nwread_min, LeftLeadParam.nwread_max);

  if ~and(isfield(LeftLeadParam,'nwdump_min'),isfield(LeftLeadParam,'nwdump_max'))
      LeftLeadParam.nwdump_min = 0;
      LeftLeadParam.nwdump_max = 0;
  end
  fprintf(1,'\t  To dump left-lead GF data for frequency points: %4d to %4d \n', LeftLeadParam.nwdump_min, LeftLeadParam.nwdump_max);  

  if ~and(isfield(RightLeadParam,'nwread_min'),isfield(RightLeadParam,'nwread_max'))
      RightLeadParam.nwread_min = 0;
      RightLeadParam.nwread_max = 0;
  end
  fprintf(1,'\t  To read right-lead GF data for frequency points: %4d to %4d \n', RightLeadParam.nwread_min, RightLeadParam.nwread_max);

  if ~and(isfield(RightLeadParam,'nwdump_min'),isfield(RightLeadParam,'nwdump_max'))
      RightLeadParam.nwdump_min = 0;
      RightLeadParam.nwdump_max = 0;
  end
  fprintf(1,'\t  To dump right-lead GF data for frequency points: %4d to %4d \n', RightLeadParam.nwdump_min, RightLeadParam.nwdump_max);  
  
  % -------------------------------------------------------------------
  % Check acoustic sum rule
  % -------------------------------------------------------------------

  err_tol = 1E-3; % relative error tolerance

  % LeftLeadParam.MatKC = LeftLeadParam.MatKC - ...
  %     diag(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC LeftLeadParam.MatKR],2));
  % err_check = sum(abs(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC LeftLeadParam.MatKR],2)));
  err_check = norm(abs(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC LeftLeadParam.MatKR],2))./abs(diag(LeftLeadParam.MatKC)));
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant acoustic sum rule violation (Left lead bulk), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Acoustic Sum Rule check passed (Left lead bulk) \n');      
  end

  % RightLeadParam.MatKC = RightLeadParam.MatKC - ...
  %     diag(sum([RightLeadParam.MatKL RightLeadParam.MatKC RightLeadParam.MatKR],2));
  % err_check = sum(abs(sum([RightLeadParam.MatKL RightLeadParam.MatKC RightLeadParam.MatKR],2)));
  err_check = norm(abs(sum([RightLeadParam.MatKL RightLeadParam.MatKC RightLeadParam.MatKR],2))./abs(diag(RightLeadParam.MatKC)));
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant acoustic sum rule violation (Right lead bulk), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Acoustic Sum Rule check passed (Right lead bulk) \n');      
  end

  % err_check = sum(abs(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC CenterParam.MatKLC],2)));
  err_check = norm(abs(sum([LeftLeadParam.MatKL LeftLeadParam.MatKC CenterParam.MatKLC],2))./abs(diag(LeftLeadParam.MatKC)));
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant acoustic sum rule violation (Left lead to channel), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Acoustic Sum Rule check passed (Left lead to channel) \n');      
  end

  % err_check = sum(abs(sum([CenterParam.MatKRC RightLeadParam.MatKC RightLeadParam.MatKR],2)));
  err_check = norm(abs(sum([CenterParam.MatKRC RightLeadParam.MatKC RightLeadParam.MatKR],2))./abs(diag(RightLeadParam.MatKC)));
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant acoustic sum rule violation. (Right lead to channel), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Acoustic Sum Rule check passed (Right lead to channel) \n');      
  end

  for nlayer = 1:length(CenterParam.Lyr)
      % CenterParam.Lyr(nlayer).MatKC = CenterParam.Lyr(nlayer).MatKC - ...
      %     diag(sum([CenterParam.Lyr(nlayer).MatKL CenterParam.Lyr(nlayer).MatKC CenterParam.Lyr(nlayer).MatKR],2));    
      % err_check = sum(abs(sum([CenterParam.Lyr(nlayer).MatKL CenterParam.Lyr(nlayer).MatKC CenterParam.Lyr(nlayer).MatKR],2)));
      err_check = norm(abs(sum([CenterParam.Lyr(nlayer).MatKL CenterParam.Lyr(nlayer).MatKC CenterParam.Lyr(nlayer).MatKR],2)) ...
                  ./abs(diag(CenterParam.Lyr(nlayer).MatKC)));        
      if gt(err_check,err_tol)
          warnmsg = sprintf('++ Possible significant acoustic sum rule violation (Channel layer %d), Epsilon = %8.4e',nlayer,err_check);
          warning(warnmsg);
      else
          fprintf(1,'\t  ... Acoustic Sum Rule check passed (Channel layer %d) \n',nlayer);      
      end    
  end

  fprintf(1,'\t  Acoustic Sum Rule check done. \n');



  % -------------------------------------------------------------------
  % Check matrix symmetry
  % -------------------------------------------------------------------

  % err_check = sum(sum(abs(LeftLeadParam.MatKL - LeftLeadParam.MatKR'),2));
  err_check = norm(LeftLeadParam.MatKL - LeftLeadParam.MatKR',1)/max([norm(LeftLeadParam.MatKL,1) norm(LeftLeadParam.MatKR',1)]);
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant matrix symmetry violation (Left lead bulk), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Matrix Symmetry check passed (Left lead bulk KL and KR) \n');      
  end

  % err_check = sum(sum(abs(LeftLeadParam.MatKC - LeftLeadParam.MatKC'),2));
  err_check = norm(LeftLeadParam.MatKC - LeftLeadParam.MatKC',1)/norm(LeftLeadParam.MatKC,1);
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant matrix symmetry violation (Left lead bulk), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Matrix Symmetry check passed (Left lead bulk KC) \n');      
  end

  % err_check = sum(sum(abs(RightLeadParam.MatKL - RightLeadParam.MatKR'),2));
  err_check = norm(RightLeadParam.MatKL - RightLeadParam.MatKR',1)/max([norm(RightLeadParam.MatKL,1) norm(RightLeadParam.MatKR',1)]);
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant matrix symmetry violation (Right lead bulk), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Matrix Symmetry check passed (Right lead bulk KL and KR) \n');            
  end

  % err_check = sum(sum(abs(RightLeadParam.MatKC - RightLeadParam.MatKC'),2));
  err_check = norm(RightLeadParam.MatKC - RightLeadParam.MatKC',1)/norm(RightLeadParam.MatKC,1);
  if gt(err_check,err_tol)
      warnmsg = sprintf('++ Possible significant matrix symmetry violation (Right lead bulk), Epsilon = %8.4e',err_check);
      warning(warnmsg);
  else
      fprintf(1,'\t  ... Matrix Symmetry check passed (Right lead bulk KC) \n');            
  end

  for nlayer = 2:(numel(CenterParam.Lyr)-1)
      % err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKC-CenterParam.Lyr(nlayer).MatKC'),2));    
      err_check = norm(CenterParam.Lyr(nlayer).MatKC-CenterParam.Lyr(nlayer).MatKC',1)/norm(CenterParam.Lyr(nlayer).MatKC,1);    
      if gt(err_check,err_tol)
          warnmsg = sprintf('++ Possible significant matrix symmetry violation (Channel layer %d), Epsilon = %8.4e',nlayer,err_check);
          warning(warnmsg);
      else
          fprintf(1,'\t  ... Matrix Symmetry check passed (Channel layer %d KC) \n',nlayer);            
      end    

      % err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKL-CenterParam.Lyr(nlayer-1).MatKR'),2))
      err_check = norm(CenterParam.Lyr(nlayer).MatKL-CenterParam.Lyr(nlayer-1).MatKR',1) ... 
                  /max([norm(CenterParam.Lyr(nlayer).MatKL,1) norm(CenterParam.Lyr(nlayer-1).MatKR',1)]);        
      if gt(err_check,err_tol)
          warnmsg = sprintf('++ Possible significant matrix symmetry violation (Channel layer %d), Epsilon = %8.4e',nlayer,err_check);
          warning(warnmsg);
      else
          fprintf(1,'\t  ... Matrix Symmetry check passed (Channel layer %d KL) \n',nlayer);            
      end    
    
      % err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKR-CenterParam.Lyr(nlayer+1).MatKL'),2));
      err_check = norm(CenterParam.Lyr(nlayer).MatKR-CenterParam.Lyr(nlayer+1).MatKL',1) ...
                  /max([norm(CenterParam.Lyr(nlayer).MatKR,1) norm(CenterParam.Lyr(nlayer+1).MatKL',1)]);        
      if gt(err_check,err_tol)
          warnmsg = sprintf('++ Possible significant matrix symmetry violation (Channel layer %d), Epsilon = %8.4e',nlayer,err_check);
          warning(warnmsg);
      else
          fprintf(1,'\t  ... Matrix Symmetry check passed (Channel layer %d KR) \n',nlayer);            
      end        
  end

  for nlayer = [1 length(CenterParam.Lyr)]
      % err_check = sum(sum(abs(CenterParam.Lyr(nlayer).MatKC-CenterParam.Lyr(nlayer).MatKC'),2));
      err_check = norm(CenterParam.Lyr(nlayer).MatKC-CenterParam.Lyr(nlayer).MatKC',1)/norm(CenterParam.Lyr(nlayer).MatKC,1);        
      if gt(err_check,err_tol)
          warnmsg = sprintf('++ Possible significant matrix symmetry violation (Channel layer %d), Epsilon = %8.4e',nlayer,err_check);
          warning(warnmsg);
      else
          fprintf(1,'\t  ... Matrix Symmetry check passed (Channel layer %d KC) \n',nlayer);            
      end    
  end

  fprintf(1,'\t  Matrix Symmetry check done. \n');

  % -------------------------------------------------------------------
  % Check size of layers in leads and channel
  % -------------------------------------------------------------------


  % -------------------------------------------------------------------
  % Check inter-layer coupling matrices between leads and channel
  % -------------------------------------------------------------------

  % if and(eq(size(LeftLeadParam.MatKC,1),size(CenterParam.Lyr(1).MatKC,1)),eq(size(LeftLeadParam.MatKC,2),size(CenterParam.Lyr(1).MatKC,2)))
  if eq(sum(abs(size(LeftLeadParam.MatKC) - size(CenterParam.Lyr(1).MatKC))),0)
      err_check = sum(sum(abs(LeftLeadParam.MatKL - CenterParam.Lyr(1).MatKL),2));
      if gt(err_check,0)
          warning('++ Possible significant inter-layer coupling matrix violation (Left lead to left contact).');
      end
  else
      warning('++ Left lead and contact matrix size mismatch.');     
  end

  if eq(sum(abs(size(RightLeadParam.MatKC) - size(CenterParam.Lyr(length(Lyr)).MatKC))),0)
      err_check = sum(sum(abs(RightLeadParam.MatKR - CenterParam.Lyr(length(Lyr)).MatKR),2));
      if gt(err_check,0)
          warning('++ Possible significant inter-layer coupling matrix violation (Right lead to right contact).');
      end
  else
      warning('++ Right lead and contact matrix size mismatch.');         
  end

  %{
  K_temp = CenterParam.MatKLC( 1:length(LeftLeadParam.MatKC), 1:length(LeftLeadParam.MatKC) );
  err_check = sum(sum(abs(LeftLeadParam.MatKR - K_temp),2));
  if gt(err_check,0)
      % fprintf(1,'\t ++ Possible significant inter-layer coupling matrix violation (Left lead to channel).\n');
      warning('++ Possible significant inter-layer coupling matrix violation (Left lead to channel).');     
  end

  n2 = size(CenterParam.MatKRC,2);
  K_temp = CenterParam.MatKRC( 1:length(RightLeadParam.MatKC), (1:length(RightLeadParam.MatKC))+n2-length(RightLeadParam.MatKC) );
  err_check = sum(sum(abs(RightLeadParam.MatKL - K_temp),2));
  if gt(err_check,0)
      % fprintf(1,'\t ++ Possible significant inter-layer coupling matrix violation (Right lead to channel).\n');
      warning('++ Possible significant inter-layer coupling matrix violation (Right lead to channel).');
  end
  %}

  fprintf(1,'\t  Coupling matrices check done. \n');

  cd(CurrDir); % return to directory containing AGF code .m files
end
