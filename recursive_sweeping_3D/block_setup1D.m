function[P]=block_setup1D(h,s3,sz,MAT)
N1=sz(1);N2=sz(2);N3=sz(3);
invA=cell(1,N3);
blockL=1/(h*h)*(s3(4:2:2*N3).*s3(3:2:2*N3-1));
blockU=1/(h*h)*(s3(2:2:2*(N3-1)).*s3(3:2:2*N3-1));

cur=1:N1*N2;
invA{1}=inv(full(MAT(cur,cur)));

for n3=2:N3
    cur=N1*N2*(n3-1)+1:N1*N2*n3;
    invA{n3}=inv(MAT(cur,cur)-(blockL(n3-1)*blockU(n3-1))*invA{n3-1});
end
P={invA,blockL,blockU};
end