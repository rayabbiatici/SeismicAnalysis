function [a,v,u,fsh]=NewmarkIntegratorNL(gamma,beta,M,C,K,P,dt,Fy,Ksh);
                                           
% Initializations    
  TimeSteps = length(P);    % for Newmark
  Time = dt*[1:TimeSteps];
  
  u = zeros(TimeSteps,1);
  v = zeros(TimeSteps,1);
  a = zeros(TimeSteps,1);
  fs = zeros(TimeSteps,1);
  fsh=zeros(TimeSteps,2);
  
  tol  = 1.0e-03;           % for Newton
  maxiters = 30;
  
% Determine initial acceleration and various time-invariant parameters
  fsa=[K Ksh Fy];
  [fs(1,:),dfs(1,:),fsh(1,:)] = ClassicalPlasticity(u(1),fsa,0,fsh(1,:));  

  a(1) = inv(M) * (P(1) - C * v(1) - fs(1,:));
  A   = (1/(beta*dt)) * M + (gamma/beta) * C;
  B   = 1/(2*beta) * M + dt * (gamma/(2*beta)-1) * C;
  Kdt = A / dt;
  
% Loop over load steps
  for i=2:TimeSteps
     du = 0; % initial guess
     dP   = P(i) - P(i-1);
     Peff =  dP + A*v(i-1) + B*a(i-1);
     
     for iters =1:maxiters 
        
        % Evaluate current iterative total displacement
          u(i) = u(i-1) + du;
        
        % Evaluate the nonlinear spring
          [fs(i,:),dfs(i,:),fsh(i,:)] = ClassicalPlasticity(u(i),fsa,fs(i-1,:),fsh(i-1,:));  

        % Compute the residual 
          res = Kdt*du + fs(i,:) - fs(i-1,:) - Peff;
          
        % Check convergence
          if abs(res) < tol
             break
          elseif (iters==maxiters & abs(res)>tol)
             display(strcat('>> Residual is = ',num2str(res)));
             warning('>> Convergence failed at'); i
          end	
  
          % Compute tangent stiffness and increment
          kt = Kdt + dfs(i,:);
          du = du - inv(kt)*res;
     end  
        
     dv = (gamma/(beta*dt))*du - (gamma/beta)*v(i-1) + ...
           dt*(1-gamma/(2*beta))*a(i-1);
     da = (1/(beta*dt^2))*du - (1/(beta*dt))*v(i-1) - ...
          (1/(2*beta))*a(i-1);
     
     v(i) = v(i-1) + dv;
     a(i) = a(i-1) + da;

  end
  
  
