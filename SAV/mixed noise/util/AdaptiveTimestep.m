function [ phi,v,tau_n,Energy_new,afa,du_x,du_y,grad_u,inter_data1,inter_data2,const_tmp] = AdaptiveTimestep( par,afa,du_x,du_y,grad_u,phi,v,inter_data1,inter_data2,DtD,tau_n,phi_old,rn,Energy,f,const_tmp)
   [phi]=new_order1(par,afa,du_x,du_y,grad_u,phi,inter_data1,inter_data2,DtD,tau_n,rn,v,const_tmp);
%   [v]= new_order2(par,afa,grad_u,phi,inter_data1,inter_data2,DtD,tau_n,rn,v,f,const_tmp);
   %step1: compute u^{n+1}
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   data_term1 = (phi-v.*log((phi+par.ep)./(v+par.ep) ) );
 %  data_term = (phi-f.*log(phi+par.ep));
   inter_data1=sum(sum(data_term1));
   
   data_term2 = (f-v).^2;
   inter_data2 = sum(sum(data_term2));
   
   afa=1+par.b*(abs(DisCuNew(phi)));
   du_x=dxf(phi);
   du_y=dyf(phi);
   grad_u=sqrt(du_x.^2+du_y.^2+sqrt(par.ep));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   E1=sum(sum(afa.*grad_u))+par.lambda1/2*inter_data1++par.lambda2/2*inter_data2;   
if E1<0
    data_term3=(phi-f).^2;
    const_tmp=-par.lambda1/2*inter_data1+par.lambda1/2*sum(sum(data_term3)) +par.const;
else
    data_term3=(phi-f).^2;
    const_tmp=-par.lambda1/2*inter_data1+par.lambda1/2*sum(sum(data_term3)) +par.const;
%     const_tmp=1;
end
%const_tmp=1;
   rn= sqrt(E1 +const_tmp);
   Energy_new=sum(sum(0.5*par.ep*grad_u.^2))+E1+const_tmp;   
%     Energy_new=sum(sum(0.5*par.ep*grad_u.^2))+sum(sum(afa.*grad_u))+par.lambda/2*sum(sum((phi-f).^2));
   err =norm(Energy-Energy_new)/norm( Energy_new);
    if err<par.tol
       tau_n = max(par.tau_min,min(par.tau_max,par.rho*sqrt(par.tol/err)*tau_n));   
    else
        Energy=Energy_new;
        tau_n= max(par.tau_min,min(par.tau_max,par.rho*sqrt(par.tol/err)*tau_n));
       [ phi,v,tau_n,Energy_new,afa,du_x,du_y,grad_u,inter_data1,inter_data2,const_tmp] = AdaptiveTimestep( par,afa,du_x,du_y,grad_u,phi,v,inter_data1,inter_data2,DtD,tau_n,phi_old,rn,Energy,f,const_tmp);
    end
end

