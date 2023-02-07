function [kvec, wkvec, vkvec] = GetPhononDispersion(Param)
% Calculates the phonon dispersion (k,w) given the matrices
    
HC = Param.MatHC;
HL = Param.MatHL;
HR = Param.MatHR;
Id = eye(size(HC)); % identity matrix 
nbmax = size(HC,1); % number of phonon branches 

a_long = Param.a_long;
kvec = (0.01:0.01:2.00)*pi/a_long; % range of wave vectors 
nkmax = numel(kvec);
dk = abs(kvec(2)-kvec(1));


ww0 = 0.1*max(abs(eig(HC)));

wkvec = [];
vkvec = [];

for nk = 1:1:nkmax
    k = kvec(nk);
    Hk = exp(-1i*k*a_long)*HL + HC + exp(1i*k*a_long)*HR;
    Vk_long = -1i*a_long*exp(-1i*k*a_long)*HL + 1i*a_long*exp(1i*k*a_long)*HR;

    [Uk,wwk] = eig(Hk+ww0*Id);
    wwk = diag(wwk - ww0*Id);
    wk = sqrt(real(wwk));

    for nb = 1:1:nbmax
        vk(nb,1) = Uk(:,nb)'*Vk_long*Uk(:,nb)/(2.0*wk(nb));
    end

    [sk, ik] = sort(wk);
    vk = vk(ik);
    wk = wk(ik);
    
    wkvec = [wkvec wk];
    vkvec = [vkvec vk];
end

wkvec = real(wkvec);
vkvec = real(vkvec);
kvec = mod(kvec+pi/a_long,2*pi/a_long)-pi/a_long;

kvec = circshift(kvec,[0, 101]);
vkvec = circshift(vkvec,[0, 101]);
wkvec = circshift(wkvec,[0, 101]);

% save('TempData.mat','kvec','wkvec','vkvec','a_long');

end
