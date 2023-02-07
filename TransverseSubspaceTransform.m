function SubInputParam = TransverseSubspaceTransform(InputParam)
% NOTE: Create the matrices associated with each Fourier component in lead (left or right)

  eta_gvec = 0E-8;

  HC = InputParam.MatHC; % square matrix
  HL = InputParam.MatHL; % square matrix
  HR = InputParam.MatHR; % square matrix

  TransCell = InputParam.TransCell; % data structure for transverse cells of lead
  n_tran = numel(TransCell);        % number of transverse cells or total number of Fourier components
  n_tran1 = InputParam.n_tran1;     % number of Fourier components (transverse 1) of lead
  n_tran2 = InputParam.n_tran2;     % number of Fourier components (transverse 2) of lead
  gvec1 = (1+1i*eta_gvec)*InputParam.gvec_tran1;    % reciprocal lattice vector (transverse 1) of lead
  gvec2 = (1+1i*eta_gvec)*InputParam.gvec_tran2;    % reciprocal lattice vector (transverse 2) of lead


  % -------------------------------------------------------------------------------------------
  % 
  % -------------------------------------------------------------------------------------------
  for nt1 = 1:1:n_tran1 % loop over Fourier components (transverse 1)
      for nt2 = 1:1:n_tran2 % loop over Fourier components (transverse 2)
          kvec = double(nt1)/double(n_tran1)*gvec1 + double(nt2)/double(n_tran2)*gvec2; 
            % corresponding transverse wave vector
                 
          nrange0 = TransCell(1).nrange; % matrix indices for 'origin' subcell in transverse direction
          rvec0 = TransCell(1).rvec;     % lattice coord. for 'origin' subcell in transverse direction
          SubHC = HC(nrange0,nrange0);
          SubHL = HL(nrange0,nrange0);
          SubHR = HR(nrange0,nrange0);
          for n = 2:n_tran % loop over transverse cells 
              nrange = TransCell(n).nrange; % matrix indices for 'current' subcell in transverse direction
              rvec = TransCell(n).rvec; % lattice coord. for 'current' subcell in transverse direction
              drvec = rvec - rvec0; % lattice coord. displacement of 'current' subcell relative to 'origin' subcell    
              phasefactor = 2.0*pi*kvec*drvec'; % positive sign here becos we are doing a Fourier transform
              SubHC = SubHC + HC(nrange0,nrange)*exp(1i*phasefactor);
              SubHL = SubHL + HL(nrange0,nrange)*exp(1i*phasefactor);
              SubHR = SubHR + HR(nrange0,nrange)*exp(1i*phasefactor);
          end
          SubInputParam(nt1,nt2).MatHC = SubHC;
          SubInputParam(nt1,nt2).MatHL = SubHL;
          SubInputParam(nt1,nt2).MatHR = SubHR;
          SubInputParam(nt1,nt2).a_long = InputParam.a_long;
          SubInputParam(nt1,nt2).n_tran1 = 1;
          SubInputParam(nt1,nt2).n_tran2 = 1;
          % disp([size(SubHC) size(SubHL) size(SubHR)]);
      end
  end

end




% === DEBUGGING CODE ===
%{
TempHC = zeros(size(HC));
TempHL = zeros(size(HL));
TempHR = zeros(size(HR));

for n1 = 1:1:n_tran
    for n2 = 1:1:n_tran
        nrange1 = TransCell(n1).nrange;
        nrange2 = TransCell(n2).nrange;

        rvec1 = TransCell(n1).rvec;
        rvec2 = TransCell(n2).rvec;
        drvec = rvec1 - rvec2;

        for nt1 = 1:1:n_tran1
            for nt2 = 1:1:n_tran2
                kvec = double(nt1)/double(n_tran1)*gvec1 + double(nt2)/double(n_tran2)*gvec2;
                phasefactor = -2.0*pi*kvec*drvec';
                TempHC(nrange1,nrange2) = TempHC(nrange1,nrange2) ...
                                          + SubInputParam(nt1,nt2).MatHC * exp(1i*phasefactor);
                TempHL(nrange1,nrange2) = TempHL(nrange1,nrange2) ...
                                          + SubInputParam(nt1,nt2).MatHL * exp(1i*phasefactor);
                TempHR(nrange1,nrange2) = TempHR(nrange1,nrange2) ...
                                          + SubInputParam(nt1,nt2).MatHR * exp(1i*phasefactor);
            end
        end
    end
end
%}

% TempHC = TempHC * 1.0/double(n_tran); % DEBUG LINE
% TempHL = TempHL * 1.0/double(n_tran); % DEBUG LINE
% TempHR = TempHR * 1.0/double(n_tran); % DEBUG LINE
% save('TempData.mat','TempHC','HC','TempHL','HL','TempHR','HR','TransCell','SubInputParam','gvec1','gvec2');
% SubInputParam = 0; % DEBUG LINE


