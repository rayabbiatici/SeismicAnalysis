function c=column(data)
%converts standard csmip ground motion data (format 8f*.*) to column
%vector (1f*.*)
[n,p]=size(data);
noelem=n*p;
data=data';
c=reshape(data,noelem,1);