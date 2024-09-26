function [EffG_LR, EffG_RR, EffG_RL, EffG_LL, P] = GreensFunctionSolver(ww,Lyr)
% NOTE: The input argument should be in the form of 
%   Lyr(n).MatHL
%   Lyr(n).MatHC
%   Lyr(n).MatHR

  nmax = length(Lyr);

  [Lyr(1:nmax).MatLeftSurfG]  = deal(Lyr(1:nmax).MatHC);
  [Lyr(1:nmax).MatRightSurfG] = deal(Lyr(1:nmax).MatHC);
  [P(1:nmax).MatG]            = deal(Lyr(1:nmax).MatHC);

  if gt(nmax,1) % i.e. we have more than 1 layer
      for n = 1:nmax
          if eq(n,1)
              Lyr(n).MatLeftSurfG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC );
          else
              Lyr(n).MatLeftSurfG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC ...
                  - Lyr(n).MatHL * Lyr(n-1).MatLeftSurfG * Lyr(n-1).MatHR );
          end
      end
      
      for n = nmax:-1:1
          if eq(n,nmax)
              Lyr(n).MatRightSurfG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC );
          else
              Lyr(n).MatRightSurfG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC ...
                  - Lyr(n).MatHR * Lyr(n+1).MatRightSurfG * Lyr(n+1).MatHL );
          end
      end
      
      
      for n = 1:nmax
          if eq(n,1)            
              P(n).MatG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC ...
                  - Lyr(n).MatHR * Lyr(n+1).MatRightSurfG * Lyr(n+1).MatHL );
          elseif eq(n,nmax)            
              P(n).MatG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC ...
                  - Lyr(n).MatHL * Lyr(n-1).MatLeftSurfG * Lyr(n-1).MatHR );
          else            
              P(n).MatG = inv( ww*eye(size(Lyr(n).MatHC)) - Lyr(n).MatHC ....
                  - Lyr(n).MatHR * Lyr(n+1).MatRightSurfG * Lyr(n+1).MatHL ...
                  - Lyr(n).MatHL * Lyr(n-1).MatLeftSurfG * Lyr(n-1).MatHR );
          end
      end    
    
      %{
      OldG = P(nmax).MatG;    
      for n = (nmax-1):-1:1
          NewG = Lyr(n).MatLeftSurfG * Lyr(n).MatHR * OldG;
          OldG = NewG;
      end        
      EffG_LR = NewG;
      %}

      % === compute left edge/right edge Green's function (off-diagonal block submatrices)
      OldG = P(nmax).MatG;
      for n = (nmax-1):-1:1
          NewG = Lyr(n).MatLeftSurfG * Lyr(n).MatHR * OldG;
          OldG = NewG;
      end        
      EffG_LR = NewG;
  
      % === compute right edge/left edge Green's function (off-diagonal block submatrices)
      OldG = P(1).MatG;
      for n = 2:1:nmax
          NewG = Lyr(n).MatRightSurfG * Lyr(n).MatHL * OldG;
          OldG = NewG;
      end        
      EffG_RL = NewG;
  
      EffG_LL = P(1).MatG;
      EffG_RR = P(nmax).MatG;
  
  else % we have only 1 layer and just do the inversion    
      EffG_LR = inv( ww*eye(size(Lyr(1).MatHC)) - Lyr(1).MatHC );
      P(1).MatG = EffG_LR;
      EffG_RL = EffG_LR;
      EffG_LL = EffG_LR;
  end
end

