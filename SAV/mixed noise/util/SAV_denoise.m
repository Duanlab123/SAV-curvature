function [HR,error,Energy_iter] = SAV_denoise(f)
% parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%model parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.ep=1e-8;
par.b=0.1;
par.lambda1=0.08;
par.lambda2=2;
par.const=0;
%%%schme parameters%%%%%%%%%%%%%%%%%%%
par.tau_min=1e-3;
par.tau_max=0.2;
par.rho=0.6;
par.tol=0.006;
Time=3;
relative_tol=1e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=f;
phi=f;
afa=1+par.b*(abs(DisCuNew(phi)));
du_x=dxf(phi);
du_y=dyf(phi);
grad_u=sqrt(du_x.^2+du_y.^2+sqrt(par.ep));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_term1 = (phi-v.*log((phi+par.ep)./(v+par.ep) ) );
inter_data1=sum(sum(data_term1));
data_term2 = (f-v).^2;
inter_data2 = sum(sum(data_term2));
E1 = sum(sum(afa.*grad_u))+par.lambda1/2*inter_data1+par.lambda2/2*inter_data2;
if E1<0
    data_term3=(phi-f).^2;
    const_tmp=-par.lambda1/2*inter_data1+par.lambda1/2*sum(sum(data_term3))+par.const;
else
    data_term3=(phi-f).^2;
    const_tmp=-par.lambda1/2*inter_data1+par.lambda1/2*sum(sum(data_term3))+par.const;
 %   const_tmp=1;
end
%  const_tmp=1;
rn=sqrt(E1+const_tmp);
Energy=sum(sum(0.5*par.ep*grad_u.^2))+E1+const_tmp;

t=2*par.tau_min;
Energy_iter=zeros(round(Time/par.tau_min)+1,1);
error=zeros(round(Time/par.tau_min)+1,1);
iter=1;
Energy_iter(iter,1)=Energy;
[l1,l2]=size(phi);
DtD = abs(psf2otf([1,-1],[l1, l2])).^2 + abs(psf2otf([1;-1],[l1, l2])).^2;
tau_n=par.tau_max;

while t<Time
    phi_old=phi;     
    [phi,v,tau_n,Energy_new,afa,du_x,du_y,grad_u,inter_data1,inter_data2,const_tmp]=AdaptiveTimestep(par,afa,du_x,du_y,grad_u,phi,v,inter_data1,inter_data2,DtD,tau_n,phi_old,rn,Energy_iter(iter,1),f,const_tmp);    
    iter=iter+1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    Energy_iter(iter,1) =Energy_new;
    error(iter)=norm( Energy_iter(iter,1)- Energy_iter(iter-1,1))/norm( Energy_iter(iter-1,1));

    if error(iter) < relative_tol
        break;
    end
    t=t+tau_n;
 end
HR=phi;
end


