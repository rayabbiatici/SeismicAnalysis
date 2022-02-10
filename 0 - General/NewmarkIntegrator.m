  function [a,v,u]=NewmarkIntegrator(gamma,beta,m,c,k,p,dt)

  n=length(p);
  u=zeros(n,1);
  v=zeros(n,1);
  a=zeros(n,1);
  
  a(1)=p(1)/m;
  K = k + (gamma/(beta*dt))*c + (1/(beta*dt^2))*m;
  A = (1/(beta*dt))*m + (gamma/beta)*c;
  B = (1/(2*beta))*m + dt*(gamma/(2*beta)-1)*c;
  
  for i=2:n
     dp = p(i)-p(i-1) + A*v(i-1) + B*a(i-1);
     du = dp/K;
     dv = (gamma/(beta*dt))*du - (gamma/beta)*v(i-1) + ...
         dt*(1-gamma/(2*beta))*a(i-1);
     da = (1/(beta*dt^2))*du - (1/(beta*dt))*v(i-1) - ...
         (1/(2*beta))*a(i-1);
     
     u(i)=u(i-1)+du;
     v(i)=v(i-1)+dv;
     a(i)=a(i-1)+da;
  end

  