%NuclearNrMinNoisy rankone
% matrix completion
close all
clearvars
Interval=20:20:50;

for rate=Interval
    rate
    for k=1:10
        k
        % Generating read matrix
        l=30;N=30;
        h=2*(rand(l,1)>.5)-1;
        v=2*(rand(l,1)>.5)-1;
        M=v*h';
        Omeg= randperm(N*l,round(rate*N*l/100)); %    randi([1 N*l],1,round(rate*N*l/100)) is not unique
        R=zeros(size(M));
        for i=1:length(Omeg)
            R(Omeg(i))=M(Omeg(i));
        end
        
        

      % Generating quality  E and error matrix Z
        E=zeros(size(R));
        lambda=4;  % a metric for mizane noise   2> 126  4>28
        aa=poissrnd(lambda,length(Omeg),1);
        
        E(Omeg)=10.^(-aa);
        A=zeros(size(R));
        Z=zeros(size(R));
        for i=1:length(Omeg)
        A(Omeg(i)) = binornd(1, E(Omeg(i)) );
        if A(Omeg(i))==1
            Z(Omeg(i))=-2*M(Omeg(i));
        end
        end
        R_ex=R;
        R=R+Z;
        
        
        
        
       
        % R is the measurment matrix and M is the unkown matrix, E is given
        omg=find(R);
        %No weight
        W=zeros(size(R))+1;
        cvx_begin quiet
        variable X(N,l)
        minimize norm_nuc(X)
        subject to
        W(omg)'*(X(omg)-R(omg)).^2 <= .1; %norm(X(omg)-R(omg),2)<= .001;
        cvx_end
        M_ht=full(X);
        e_Nuc(k)=norm(M_ht-M);
        
       
        [U sig V]=svd(M_ht);
        M_ht_svd= U(:,1)*sig(1,1)*V(:,1)';
        e_Nuc_svd(k)=norm(M_ht_svd-M);
        M_ht_svd_round=2*(M_ht_svd>0)-1;
        e_Nuc_svd_round(k)=norm(M_ht_svd_round-M);


        
        
        
        %weighted
        W=zeros(size(R));
        W(omg)=log2(1./E(omg));
        cvx_begin quiet
        variable X(N,l)
        minimize norm_nuc(X)
        subject to
        W(omg)'*(X(omg)-R(omg)).^2 <= .1; %norm(X(omg)-R(omg),2)<= .001;
        cvx_end
        Mw_ht=full(X);
        ew_Nuc(k)=norm(Mw_ht-M);
      
        
        
        [Uw sigw Vw]=svd(Mw_ht);
        Mw_ht_svd= Uw(:,1)*sigw(1,1)*Vw(:,1)';
        ew_Nuc_svd(k)=norm(Mw_ht_svd-M);
        Mw_ht_svd_round=2*(Mw_ht_svd>0)-1;
        ew_Nuc_svd_round(k)=norm(Mw_ht_svd_round-M);

        
        
        
        
    end
    
e_Nuc_m(rate)=mean(e_Nuc);
        e_Nuc_svd_m(rate)=mean(e_Nuc_svd);
        e_Nuc_svd_round_m(rate)=mean(e_Nuc_svd_round);
    
        
ew_Nuc_m(rate)=mean(ew_Nuc);
        ew_Nuc_svd_m(rate)=mean(ew_Nuc_svd);
        ew_Nuc_svd_round_m(rate)=mean(ew_Nuc_svd_round);
    
    
end



%10*log10(mean(e))-10*log10(mean(e_w))

figure(1)
hold on
plot(Interval,(e_Nuc_svd_round_m(Interval)),Interval,(e_Nuc_svd_m(Interval)),Interval,(e_Nuc_m(Interval)),'LineWidth',2)
%legend()
title(' l=40 N=40   delta=.1')
xlabel('rate of known entries')
ylabel(' norm of error matrix')



plot(Interval,(ew_Nuc_svd_round_m(Interval)),Interval,(ew_Nuc_svd_m(Interval)),Interval,(ew_Nuc_m(Interval)),'LineWidth',2)
legend('round','svd','nuc','round -w ','svd -w ','nuc-w ')
title(' l=40 N=40   delta=.1')
xlabel('rate of known entries')
ylabel(' norm of error matrix')




