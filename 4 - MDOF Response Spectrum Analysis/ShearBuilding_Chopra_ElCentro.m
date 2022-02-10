% **********************************************************************
% *                                                                    *
% *     This program determines the response of a multistory shear     *
% *     building (infinitely rigid slabs) under a given ground motion. *
% *                                                                    *
% *                     |                 |                            *
% *                     O======( mi )=====O                            *
% *                     |                 |                            *
% *            0.5*k_i  |                 | 0.5*k_i                    *
% *                     |                 |                            *
% *                     O=====( m_i-1)====O                            *
% *                     |                 |                            *
% *                                                                    *
% *                                                                    *
% *                     |                 |                            *
% *                     O======( m1 )=====O                            *
% *                     |                 |                            *
% *            0.5*k_1  |                 | 0.5*k_1                    *
% *                     |                 |                            *
% *                   ~~~~~             ~~~~~                          *
% *               <---------->      <---------->                       *
% *                   ag(t)             ag(t)                          *
% *                                                                    *
% *-----INPUT:                                                         *
% *                                                                    *
% *      <dt>    : Time increment                                      *
% *      <ag>    : nsteps x 1 vector of ground motion accelerations    *
% *      <kfl>   : ndof x 1 vector of floor (bending) stiffnesses      *
% *      <mfl>   : ndof x 1 vector of floor masses                     *
% *      <kdamp> : Stiffness proportional damping coefficient          *
% *      <mdamp> : Mass proportional damping coefficient               *
% *                                                                    *
% *-----OUTPUT:                                                        *
% *                                                                    *
% *      <u,v,a> : ndof x nsteps matrices of floor displacements,      *
% *                velocities and accelerations                        *
% *      <M,C,K> : ndof x ndof matrices of mass, damping and stiffness *
% *      <w>     : ndof x 1 vector of modal frequencies                *
% *      <xi>    : ndof x 1 vector of modal damping ratios             *
% *      <Phi>   : ndof x ndof matrix of mode shapes                   *
% *                                                                    *
% *-----GLOSSARY:                                                      *
% *                                                                    *
% *      <ndof>  : Number of degrees-of-freedom (floors)               *
% *                                                                    *
% **********************************************************************
% *                                                                    *
% *   ShearBuilding                                                    *
% *                                                                    *
% *                                                                    *
% *                                                                    *
% **********************************************************************


% Initializations
  clf
  clc
  clear all
 
% Give the project a title
  ptitle =  ' Chopra Five-Story Shear Building'

% El Centro NS ground motion
 
  load ElCentroNS.txt
  ag90 = column (ElCentroNS);
     
% Set units for ground motion
  ag = 386.4 * ag90;  % north-south acceleration                  
  clear ag90 ; clear El_Centro_Chopra_194x8.txt; % Remove raw data from memory
  dt = 1/50; % Time-Step  
%   time = dt:dt:dt*length(ag);  % Time Array
  time = 0:dt:dt*(length(ag)-1); %Time Span

  % Number of stories
  nStories = 5;
 
  % Floor heights
  h1  = 12.0; % Height
  hfl = h1*ones(1, nStories);
  Htot = nStories*h1 % Total Height of the Structure...
 
% Floor lateral stiffnesses
  EI1 = 4542.8; k1 = 12*EI1/h1^3;  % gives k=31.54 k/in
  kfl = k1*ones(1, nStories); % inserting the stiffness per floor
 
% Floor masses
  m1  = 100/386.4; % 100-kips/g...to determine mass
  mfl = m1*ones(1, nStories); % pre-allocating the length of the arrray in accordance to the structure. mass per story
  Mtot = nStories*m1 % total mass of the structure

 % Pre-allocation
  ndof = length(mfl)    % Number of degrees-of-freedom
  M    = zeros(ndof);   C    = zeros(ndof);  
  K    = zeros(ndof);   Phi  = zeros(ndof);
  w    = zeros(ndof,1); xi   = zeros(ndof,1);
  s    = zeros(ndof,ndof);
 
  D    = zeros(ndof,length(ag));
  A    = zeros(ndof,length(ag));
 
% Construct system matrices
  for i=1:ndof-1  
     % creating the diagonals for both the mass and stiffness matrix...
     K(i,i)   = kfl(i) + kfl(i+1);
     K(i+1,i) = -kfl(i); K(i,i+1) = -kfl(i);
     M(i,i)   = mfl(i);
  end
 
  % Last value in there...
  K(ndof,ndof) = kfl(ndof);
  M(ndof,ndof) = mfl(ndof);
 
% Solve eigenvalue problem for natural frequencies and mode shapes
  [Phi,wsq] = eig(K,M); %eig function...
 
% Determine modal frequencies

  w=sqrt(diag(wsq));

  for j=1:nStories
   T(j)  = (2*pi)/w(j);% converting each natural frequency into a period...
  end
% Determine modal damping ratios

% Rayleigh Damping coefficients and matrix
  dratio=0.02 ; % Set damping to 2% to allow comparison with Chopra
  a1 = 2*dratio/(w(1)+w(4)) ; % 2% @ mode 1&4
  a0 = a1*w(1)*w(4) ; % 2% @ mode 1&4
 
  for i=1:ndof
      xi(i)=a0/(2*w(i))+a1*w(i)/2 ;
  end
 
  C = a0*M+a1*K ;
  CN = Phi'*C*Phi ;
 
  fprintf('The Rayleigh Damping coefficients for all modes:\r'); xi
  fprintf('The Rayleigh Damping coefficients a0 and a1:\r'); a0, a1
  fprintf('The Rayleigh damping Matirx:\r'); C
   
  for i=1:ndof %scanning for each dof or per mode...
     
      % Modal mass, damping, and stiffeness
      xi(i) = dratio;  
      m(i) = Phi(:,i)'*M*Phi(:,i); % "diagonalizing"
      c(i) = 2*xi(i)*w(i)*m(i) ; % this is basically the damping per mode...this will be used in the newmark beta
      k(i) = Phi(:,i)'*K*Phi(:,i); % "diagonalizing"...but not really
     
      % Participation Factors & distirbution vectors
      L(i) = Phi(:,i)'*M*ones(ndof,1); %slide L6-17
      % H(i) is the numerator of h* L6-25...array of increasing height
      H(i) = Phi(:,i)'*M*cumsum(hfl)'; %cumsum is a function that cummulatively adds and places in an array.
      Gamma(i) = L(i)/m(i);% slide L6-17
     
      % Effective modal mass & height
      meff(i) = (L(i)^2)/m(i); %slide L6-24
      MMP(i) = meff(i)/Mtot %slide L6-24
      heff(i) = H(i)/L(i); %hn* L6-25
      Heff(i) = heff(i)/Htot*100 %percentage/ ratio of the effective height w.r.t the total height.
      Heffh(i) = Heff(i)/100*nStories % w.r.t. a single story height.
      s(:,i) = Gamma(i)*M*Phi(:,i); %L6 -23...also called modal mass particiaption. MASS!!!!
  end
 
  %Difference between fn, sn, and rn....Chopra 516

% The natural frequencies, mode shapes, and effective masses & heights
  fprintf('The modal masses, m, are:\r'); m
  fprintf('The modal frequencies, w, are:\r'); w
  fprintf('The modal periods, T, are:\r'); T'
  fprintf('The mode shapes, phi, are;\r'); Phi
  fprintf('The damping matrix, C, is:\r'); C
  fprintf('The modal damping ratios, xi, are:\r'); xi
  fprintf('The Modal Participation factors are:\r'); Gamma
  fprintf('The Effective modal masses, meff, are:\r'); meff
  fprintf('The Modal Mass Participation, MMP, are:\r'); MMP
  fprintf('The Effective modal heights, heff, are:\r'); heff
  fprintf('The Effective modal height(%), Heff/Htot, are:\r'); Heff
  fprintf('The Effective modal height(%), Heff/Htot, are:\r'); Heffh
   
% Compute response of the building
  NewmarkGamma = 1/2; NewmarkBeta = 1/6;     % Linear acceleration method

  for i=1:ndof
     
     P = -m(i)*ag';
     % the m, c, and k are all de-coupled...slide L6-27...L6-16 to L6-19
     [a, v, u] = NewmarkIntegrator(NewmarkGamma,NewmarkBeta,m(i),c(i),k(i),P,dt);
   
     D(i,:) = u;%what we really need.
     A(i,:) = w(i)^2*D(i,:); % This is psuedo-acceleration
   
  end
 
% Total Response
  u=Phi*diag(Gamma)*D; %L6-27...total response of the structure throughout time of the gm...
                       %also slide L6-19.... basically same equation
 
% Decomposition of top floor disp
  for i=1:ndof
      rst(:,i) = inv(K)*s(:,i);% Chopra 516... in modal terms...Nara 4-44...aka displacement vector...xn...please clarify
      udecomp = rst(:,i)*A(i,:);% since (PSA)*(1/(wn^2)) = D... = xn
      utop(i,:) = udecomp(ndof,:);% displacement of top floor...for all the modes...but only for the top floor.
  end
 
  d1tmax = max(utop,[],2) % what does this do again?
  dltmin = min(utop,[],2) % also this one?
 
% Response synthesis matrix  
  B      = zeros(4,ndof); % pre-allocation
  B(1,:) = zeros(1, ndof);       % Top story displacement
  B(1,ndof) = 1;
  B(2,:) = zeros(1, ndof);      % Story drift between 3rd and 4th floors
  B(2,3) = -1; B(2,4) = 1;
  B(3,:) = (K*ones(1, ndof)')'; % Base shear
  B(4,:) = (K*cumsum(hfl)')'; % Overturning moment
  R = B*u;
   
% Plot of the mode shapes
  color = hsv(ndof);  level = [0 1:ndof];  figure(1);
  for i = 1:ndof
      plot([0; Phi(:,i)],level,'color',color(i,:));
      hold on;
  end
  grid on
  clear leg
%   for i=1:ndof
%      leg(i,:) = strcat(['Mode ',num2str(i),' ']);
%   end  
%   legend(leg)
  title('Mass Normalized Mode Shapes')
  xlabel('Displacement');  ylabel('Level');  hold off


% Plot of response quantities, R
   figure(2);  
   subplot(4,1,1); plot(time,R(1,:),'r');
   title('Top Displacement vs. Time')
   ylabel('Displacement (in)')
   grid on
   
   subplot(4,1,2); plot(time,R(2,:),'g');
   title('Interstory Drift Between 3rd and 4th Floors vs. Time')
   ylabel('Drift (in)')
   grid on
   
   subplot(4,1,3); plot(time,R(3,:),'b');
   title('Base Shear vs. Time')
   ylabel('Base Shear (kips)')
   grid on
   
   subplot(4,1,4); plot(time,R(4,:),'r');
   title('Overturning Moment vs. Time')
   xlabel('Time (sec)'); ylabel('Moment (kip-ft)')
   grid on
   
 % Plot modal contributions to top displacement
 
   figure(3);
   subplot(5,1,1); plot(time,utop(1,:));
   title('Modal Contribution to Top Displacement')
   ylabel('1st Mode')
   axis([0 30 -10.0 10.0])
   grid on
   
   subplot(5,1,2); plot(time,utop(2,:));
   ylabel('2nd Mode')
   axis([0 30 -2.5 2.5])
   grid on
   
   subplot(5,1,3); plot(time,utop(3,:));
   ylabel('3rd Mode')
   axis([0 30 -2.5 2.5])
   grid on
   
   subplot(5,1,4); plot(time,utop(4,:));
   ylabel('4th Mode')
   axis([0 30 -2.5 2.5])
   grid on
   
   subplot(5,1,5); plot(time,utop(ndof,:));
   ylabel('5th mode')
   xlabel('Time (sec)')
   axis([0 30 -2.5 2.5])
   grid onUpd
   
   % Answer Check Plot Vb & Mb from Response Synthesis Matrix and results
   % from using effective mass & heights, should be identical.
   
   for i=1:ndof
       vb(i,:) = meff(i)*A(i,:);
       mb(i,:) = meff(i)*heff(i)*A(i,:);
   end
   
   Vb = sum(vb);
   Mb = sum(mb);
   
   figure(4)
   subplot(2,1,1)
   plot(time,Vb, 'Linewidth', 2); hold on;
   plot(time,R(3,:), 'r--', 'Linewidth', 2); hold off;
   title('Base Shear vs. Time')
   ylabel('Base Shear (kips)')
   grid on
   
   subplot(2,1,2)
   plot(time,Mb, 'Linewidth', 2); hold on;
   plot(time,R(4,:), 'r--', 'Linewidth', 2); hold off;
   title('Overturning Moment vs. Time')
   xlabel('Time (sec)'); ylabel('Moment (kip-ft)')
   grid on
