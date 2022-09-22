function [phi]=new_order1(par,afa,du_x,du_y,grad_u,phi,inter_data1,inter_data2,DtD,tau_n,rn,v,const_tmp)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   E1=sum(sum(afa.*grad_u))+par.lambda1/2*inter_data1 + par.lambda2/2*inter_data2;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   div=dxb(du_x./sqrt(grad_u.^2+par.ep))+dyb(du_y./sqrt(grad_u.^2+par.ep));
   dF=-afa.*div+(par.lambda1./(phi+par.ep)).*(phi-v);
   bn=dF./sqrt(E1+const_tmp);  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   temp=Ainverse(par,bn,DtD,tau_n);
   gamma1n=sum(sum(temp.*bn));
   %%(bn,u_n+1)
   cn=phi-tau_n*rn.*bn+0.5*tau_n*bn*sum(sum(bn.*phi));%cn
   tempg=Ainverse(par,cn,DtD,tau_n);%A^(-1)cn
   gamma2n=sum(sum(bn.*tempg));%(bn,A^(-1)cn)
   tempf=gamma2n/(1+0.5*tau_n*gamma1n);%(bn,u_n+1)
   %%u_n+1
   phi=tempg-0.5*tau_n*temp.*tempf;
end