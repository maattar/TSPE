function d_d = a2d(d_a)

dim_a = size(d_a,1);
dim_d = 1+((dim_a-1)/10);

d_d = NaN(dim_d,size(d_a,2));

for i=1:dim_d-1
    d_d(i,:) = sum(d_a(10*i-9:10*i,:))/10;
end
