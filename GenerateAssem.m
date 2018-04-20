
function [Sr,Sx,Sy,S] = GenerateAssem(rho, Ma, Mb)

% given the measurements Ma for Alice, Mb for Bob, and the state rho for
% Alice-Bob-Charlie, generates the assemblage for Charlie in CG form

    [da,~,ma] = size(Ma);
    [db,~,mb] = size(Mb);
    dtot=size(rho,1);
    dc=dtot/(da*db);
    
    Sr = PartialTrace(rho,[1,2],[da,db,dc]);
    
    Sx=zeros(dc,dc,ma);
    for x=1:ma
        Sx(:,:,x)=PartialTrace(Tensor(Ma(:,:,x),eye(db),eye(dc))*rho,[1,2],[da,db,dc]);
    end

    Sy=zeros(dc,dc,mb);
    for y=1:mb
        Sy(:,:,y)=PartialTrace(Tensor(eye(da),Mb(:,:,y),eye(dc))*rho,[1,2],[da,db,dc]);
    end
    
    S=zeros(dc,dc,ma*mb);
    for x=1:ma
        for y=1:mb
            S(:,:,ma*(y-1)+x)=PartialTrace(Tensor(Ma(:,:,x),Mb(:,:,y),eye(dc))*rho,[1,2],[da,db,dc]);
        end
    end
    
end