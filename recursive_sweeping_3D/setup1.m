function[P]=setup1(NPML,NLPD,NPAD,pL,pR,h,ksq,s1,s2,s3)
[N1,N2,N3]=size(ksq);

%generate the partitions
M=(N1+1)/2;
ordL=NLPD+NPML-1:NLPD:M-1;
NP=length(ordL);
pttnL=cell(1,NP);
pttnL{1}=[1,NLPD+NPML-1];
for g=2:NP
    pttnL{g}=[ordL(g-1)+1,ordL(g)];
end
ordR=N1+1-ordL(end:-1:1);
pttnR=cell(1,NP);
for g=1:NP-1
    pttnR{g}=[ordR(g),ordR(g+1)-1];
end
pttnR{NP}=[N1+1-(NLPD+NPML-1),N1];
pttnM=[ordL(end)+1,ordR(1)-1];

    %construct left sweeping matrices
    PL=cell(1,NP);
    if(1)
        b=1;
        fm=pttnL{b}(1);to=pttnL{b}(2);
        ksqnew=ksq(fm:to,:,:);
        s1new=s1(2*fm-1:2*to+1);
        PL{b}=setup2(NPML,NLPD,NPAD,pL,pR,h,ksqnew,s1new,s2,s3);
        fprintf('complete 1: PL fm=%d to=%d\n',fm,to);
    end
    for b=2:NP
        fm=pttnL{b}(1);to=pttnL{b}(2);
        ksqnew=ksq(fm-(NPAD-1):to,:,:);
        s1new=[pL,s1(2*fm:2*to+1)];
        PL{b}=setup2(NPML,NLPD,NPAD,pL,pR,h,ksqnew,s1new,s2,s3);
        fprintf('complete 1: PL fm=%d to=%d\n',fm,to);
    end
    %construct right sweeping matrices
    PR=cell(1,NP);
    for b=1:NP-1
        fm=pttnR{b}(1);to=pttnR{b}(2);
        ksqnew=ksq(fm:to+(NPAD-1),:,:);
        s1new=[s1(2*fm-1:2*to),pR];
        PR{b}=setup2(NPML,NLPD,NPAD,pL,pR,h,ksqnew,s1new,s2,s3);
        fprintf('complete 1: PR fm=%d to=%d\n',fm,to);
    end
    if(1)
        b=NP;
        fm=pttnR{b}(1);to=pttnR{b}(2);
        ksqnew=ksq(fm:to,:,:);
        s1new=s1(2*fm-1:2*to+1);
        PR{b}=setup2(NPML,NLPD,NPAD,pL,pR,h,ksqnew,s1new,s2,s3);
        fprintf('complete 1: PR fm=%d to=%d\n',fm,to);
    end
    %construct middle layers
    if(1)
        fm=pttnM(1);to=pttnM(2);
        ksqnew=ksq(fm-(NPAD-1):to+(NPAD-1),:,:);
        s1new=[pL,s1(2*fm:2*to),pR];
        PM=setup2(NPML,NLPD,NPAD,pL,pR,h,ksqnew,s1new,s2,s3);
        fprintf('complete 1: PM fm=%d to=%d\n',fm,to);
    end

    P={{PL,pttnL},{PR,pttnR},{PM,pttnM},size(ksq)};
end