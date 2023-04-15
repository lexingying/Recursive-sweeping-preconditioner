function[u]=block_solve1D(sz,P,f)
N1=sz(1);N2=sz(2);N3=sz(3);
invA=P{1};blockL=P{2};blockU=P{3};

u=f(:);
cur=1:N1*N2;
u(cur)=invA{1}*u(cur);

for n3=2:N3
    prev=cur;
    cur=N1*N2*(n3-1)+1:N1*N2*n3;
    u(cur)=invA{n3}*(u(cur)-blockL(n3-1)*u(prev));
end

for n3=N3-1:-1:1
    prev=cur;
    cur=N1*N2*(n3-1)+1:N1*N2*n3;
    u(cur)=u(cur)-invA{n3}*(blockU(n3)*u(prev));
end
u=reshape(u,[N1,N2,N3]);
end