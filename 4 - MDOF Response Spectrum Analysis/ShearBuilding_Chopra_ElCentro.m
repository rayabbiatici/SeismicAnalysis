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
% *                                                                    *
% **********************************************************************


  
% Give the project a title
  ptitle =  ' Chopra Five-Story Shear Building'

% El Centro NS ground motion
  
  load ElCentroNS.txt 
  ag90 = column (ElCentroNS);
     
% Set units for ground motion
  ag = 386.4 * ag90;  % north-south acceleration                   
  clear ag90 ; clear ElCentroNS.txt; % Remove raw data from memory
  dt = 1/50;   
  time = dt:dt:dt*length(ag);  % Time Array
  
  % Floor heights
  h1  = 12.0;
  hfl = [h1, h1, h1, h1, h1];
  Htot = 5*h1
  
% Floor lateral stiffnesses
  EI1 = 4542.8; k1 = 12*EI1/h1^3;  % gives k=31.54 k/in
  kfl = [k1, k1, k1, k1, k1];
  
% Floor masses 
  m1  = 100/386.4;  
  mfl = [m1, m1, m1, m1, m1];
  Mtot = 5*m1

  ndof = length(mfl)    % Number of degrees-of-freedom
  M    = zeros(ndof);   C    = zeros(ndof);  
  K    = zeros(ndof);   Phi  = zeros(ndof);
  w    = zeros(ndof,1); xi   = zeros(ndof,1);
  s    = zeros(ndof,ndof);
  
  D    = zeros(ndof,length(ag));
  A    = zeros(ndof,length(ag));
  
% Construct system matrices
  for i=1:ndof-1   
     K(i,i)   = kfl(i) + kfl(i+1);
     K(i+1,i) = -kfl(i); K(i,i+1) = -kfl(i);
     M(i,i)   = mfl(i);
  end
  
  K(ndof,ndof) = kfl(ndof)
  M(ndof,ndof) = mfl(ndof)
  
% Solve eigenvalue problem for natural frequencies and mode shapes
  [Phi,wsq] = eig(K,M); 
  
% Determine modal frequencies

  w=sqrt(diag(wsq));

  for j=1:5
   T(j)  = (2*pi)/w(j);
  end
% Determine modal damping ratios

% Rayleigh Damping coefficients and matrix
  dratio=0.02 ;
  a1 = 2*dratio/(w(1)+w(4)) ;
  a0 = a1*w(1)*w(4) ;
  
  for i=1:ndof
      xi(i)=a0/(2*w(i))+a1*w(i)/2 ;
  end
 
  C = a0*M+a1*K ;
  CN = Phi'*C*Phi ;
  cn = diag(CN);
  
  fprintf('The Rayleigh Damping coefficients for all modes:\r'); xi
  fprintf('The Rayleigh Damping coefficients a0 and a1:\r'); a0, a1
  fprintf('The Rayleigh damping Matirx:\r'); C
  fprintf('The modal damping matirx PhiT*C*phi:\r'); cn
   
  for i=1:ndof
      
      % Modal mass, damping, and stiffeness
      xi(i) = 0.02;  % Set damping to 2% to allow comparison with Chopra 
      m(i) = Phi(:,i)'*M*Phi(:,i);
      c(i) = 2*xi(i)*w(i)*m(i) ;
      k(i) = Phi(:,i)'*K*Phi(:,i);
      
      % Participation Factors & distirbution vectors
      L(i) = Phi(:,i)'*M*ones(ndof,1);
      H(i) = Phi(:,i)'*M*cumsum(hfl)';
      Gamma(i) = L(i)/m(i);
      
      % Effective modal mass & height
      meff(i) = (L(i)^2)/m(i);
      MMP(i) = meff(i)/Mtot
      heff(i) = H(i)/L(i); 
      Heff(i) = heff(i)/Htot*100
      Heffh(i) = Heff(i)/100*5
      s(:,i) = Gamma(i)*M*Phi(:,i);
  end
 
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

  

  
% Decomposition of top floor disp u5(t) and max/min
% Below, utop(mode,time) is used
  
% Response synthesis matrix: top displ, story drift (3rd - 4th), Vb, Mo 
 
   
% Plot of the mode shapes
  color = hsv(ndof);  level = [0 1:ndof];  figure(1);
  for i = 1:ndof;
      plot([0; Phi(:,i)],level,'color',color(i,:));
      hold on;
  end
  grid on
  clear leg
  for i=1:ndof
     leg(i,:) = strcat('Mode ',num2str(i),' ');
  end  
  legend(leg)
  title('Mass Normalized Mode Shapes')
  xlabel('Displacement');  ylabel('Level');  hold off


% Plot of response quantities, R matirx is used
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
   
   subplot(5,1,5); plot(time,utop(5,:));
   ylabel('5th mode')
   xlabel('Time (sec)') 
   axis([0 30 -2.5 2.5])
   grid on
   
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
   plot(time,Vb,time,R(3,:))
   title('Base Shear vs. Time')
   ylabel('Base Shear (kips)')
   grid on
   
   subplot(2,1,2)
   plot(time,Mb,time,R(4,:))
   title('Overturning Moment vs. Time')
   xlabel('Time (sec)'); ylabel('Moment (kip-ft)')
   grid on
