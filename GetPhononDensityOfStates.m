function PhononDOS = GetPhononDensityOfStates(w,LeadPhonon)
% NOTE: This function computes the un-normalized phonon density of states for the lead
% NOTE: The normalization factor should be 1/Vol 
% NOTE: In 2D, Vol is the 2D area of the AGF principal layer 
% NOTE: In 3D, Vol is the 3D volume of the AGF principal layer 
  G = LeadPhonon.MatBulkG;
  PhononDOS = 1i*w*trace(G-G')/pi;
end
