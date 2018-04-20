
function [output,Fr,Fx,Fy,F,beta,Gamma] = IsAQAssemblage(Sr,Sx,Sy,S,ma,mb)

% given an assemblage in CG notation, this function determines whether it belongs to the 
% Almost-Quantum outer approximation to the set of quantum assemblages. 
% This function returns out=0 if the assemblage is outside and out=1 if inside.
%
% Steering scenario with two uncharacterised parties (Alice and Bob), and
% one characterised one (Charlie).
%
% mb = number of measurements for Alice
% mc = number of measurements for Bob
%
% this code assumes that all measurements are dichotomic


d=size(Sr,1); 
dimG=(ma+1)*(mb+1)*d; 

cvx_begin sdp quiet

    variable Gamma(dimG,dimG) hermitian % moment matrix
    
    variable pi_r(d,d) hermitian 
    variable pix(d,d,ma) hermitian 
    variable piy(d,d,mb) hermitian 
    variable pixy(d,d,ma*mb) hermitian semidefinite

    dual variable Fr
    dual variables Fx{ma}
    dual variables Fy{mb}
    dual variables F{ma*mb}
            
    minimise trace(pi_r)
            
    Gamma == hermitian_semidefinite(dimG);
                
    for x = 1:ma
        for y = 1:mb
            pix(:,:,x)-pixy(:,:,x+ma*(y-1)) == hermitian_semidefinite(d);
            piy(:,:,y)-pixy(:,:,x+ma*(y-1)) == hermitian_semidefinite(d);
            pi_r-pix(:,:,x)-piy(:,:,y)+pixy(:,:,x+ma*(y-1)) == hermitian_semidefinite(d);
        end
    end
            
    % First row
    Fr : Gamma(1:d,1:d)==Sr+pi_r;
                
    % Alice's marginal assemblage
    for k=1:ma
        Fx{k} : Gamma(1:d,k*d+1:(k+1)*d)==Sx(:,:,k)+pix(:,:,k);
    end
                
    % Bob's marginal assemblage
    for k=1:mb
        Fy{k} : Gamma(1:d,k*d*(ma+1)+1:k*d*(ma+1)+d)==Sy(:,:,k) + piy(:,:,k);
    end

    % the assemblage for all-but-one outcome
    for k=1:ma
        for j=1:mb
            F{ma*(j-1)+k} : Gamma(1:d,j*d*(ma+1)+k*d+1:j*d*(ma+1)+k*d+d)...
                == S(:,:,ma*(j-1)+k) + pixy(:,:,ma*(j-1)+k);
        end
    end
            
    % Conditions on Gamma to be AQ                
    for y=0:ma
        for z=0:mb
            for yp=0:ma
                for zp=0:mb 
                    if y==0 && z==0 % reductions of type AB, where it can be A=1 and/or B=1
                        r1=1;
                        c1=zp*d*(ma+1)+yp*d+1;
                        r2=yp*d+1;
                        c2=zp*d*(ma+1)+1;
                        r3=yp*d+1;
                        c3=zp*d*(ma+1)+yp*d+1;
                        r4=zp*d*(ma+1)+1;
                        c4=zp*d*(ma+1)+yp*d+1;
                        r5=zp*d*(ma+1)+yp*d+1;
                        c5=zp*d*(ma+1)+1;
                        r6=zp*d*(ma+1)+yp*d+1;
                        c6=zp*d*(ma+1)+yp*d+1;
                        r7=zp*d*(ma+1)+1;
                        c7=yp*d+1;
                        r8=zp*d*(ma+1)+yp*d+1;
                        c8=1;
                        r9=zp*d*(ma+1)+yp*d+1;
                        c9=yp*d+1;
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r2:r2+d-1,c2:c2+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r3:r3+d-1,c3:c3+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r4:r4+d-1,c4:c4+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r5:r5+d-1,c5:c5+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r6:r6+d-1,c6:c6+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r7:r7+d-1,c7:c7+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r8:r8+d-1,c8:c8+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r9:r9+d-1,c9:c9+d-1);
                    end
                    
                    if y==0 && z~=0 && yp~=0 && zp~=0 % Reduction of type ABB', none = 1. 
                        r1=z*d*(ma+1)+1;
                        c1=zp*d*(ma+1)+yp*d+1;
                        r2=z*d*(ma+1)+yp*d+1;
                        c2=zp*d*(ma+1)+1;
                        r3=z*d*(ma+1)+yp*d+1;
                        c3=zp*d*(ma+1)+yp*d+1;
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r2:r2+d-1,c2:c2+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r3:r3+d-1,c3:c3+d-1);
                    end
                    
                    if z~=0 && yp==0 && zp~=0 % (AB,A')=(A'B,A)', B can be =1.
                        if z~=zp
                            r1=z*d*(ma+1)+y*d+1;
                            c1=zp*d*(ma+1)+1;
                            r2=zp*d*(ma+1)+y*d+1;
                            c2=z*d*(ma+1)+1;
                            Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r2:r2+d-1,c2:c2+d-1)';
                        end
                    end
                    
                    if z==0 && y~=0 && yp~=0 && zp~=0 % Reduction of type AA'B, none = 1.
                        r1=y*d+1;
                        c1=zp*d*(ma+1)+yp*d+1;
                        r2=zp*d*(ma+1)+y*d+1;
                        c2=yp*d+1;
                        r3=zp*d*(ma+1)+y*d+1;
                        c3=zp*d*(ma+1)+yp*d+1;
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r2:r2+d-1,c2:c2+d-1);
                        Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r3:r3+d-1,c3:c3+d-1);
                    end
                    
                    if y~=0 && yp~=0 && zp==0 % (AB,B')=(AB',B)'
                        if y~=yp
                            r1=z*d*(ma+1)+y*d+1;
                            c1=yp*d+1;
                            r2=z*d*(ma+1)+yp*d+1;
                            c2=y*d+1;
                            Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r2:r2+d-1,c2:c2+d-1)';
                        end                        
                    end
                    
                    if y~=0 && z~=0 && yp~=0 && zp==0 % (AB, A'B') = (A'B',AB)'
                        if y~=yp && z~=zp
                            r1=z*d*(ma+1)+y*d+1;
                            c1=zp*d*(ma+1)+yp*d+1;
                            r2=zp*d*(ma+1)+yp*d+1;
                            c2=z*d*(ma+1)+y*d+1;
                            Gamma(r1:r1+d-1,c1:c1+d-1)==Gamma(r2:r2+d-1,c2:c2+d-1)';
                        end                                            
                    end
                end
            end
        end
    end
    
    cvx_end
    
    threshold = 1e-8; % threshold for being 'inside' AQ (i.e. anything smaller will be considered inside, anything larger outside)
    
    output=cvx_optval <= threshold;
    beta = cvx_optval;
    
end