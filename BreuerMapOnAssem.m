
function [Rr, Rx, Ry, R] = BreuerMapOnAssem(U, Sr, Sx, Sy,S,ma,mb)
% Given an antisymmetric unitary U, it applies the map
% Lambda[rho] = 1/2*(trace(rho) Id - rho - U rho' U') 
% to the assemblage to generate a new one. We are assuming a steering 
% scenario with two unstrusted parties, Alice and Bob,  ma/mb = no. 
% of measurements for Alice/Bob and the assemblage is given in CG form

    dc=size(Sr,1);
    Id=eye(dc);
    
    Rr=(trace(Sr)*Id-Sr-U*transpose(Sr)*U')/2;
    
    Rx=zeros(dc,dc,ma);
    for x=1:ma
        Rx(:,:,x)=(trace(Sx(:,:,x))*Id-Sx(:,:,x)-U*transpose(Sx(:,:,x))*U')/2;
    end
   
    Ry=zeros(dc,dc,mb);
    for y=1:mb
        Ry(:,:,y)=(trace(Sy(:,:,y))*Id-Sy(:,:,y)-U*transpose(Sy(:,:,y))*U')/2;
    end
    
    R=zeros(dc,dc,ma*mb);
    for x=1:ma
        for y=1:mb
            R(:,:,ma*(y-1)+x)=(trace(S(:,:,ma*(y-1)+x))*Id-S(:,:,ma*(y-1)+x)-U*transpose(S(:,:,ma*(y-1)+x))*U')/2;
        end
    end

end