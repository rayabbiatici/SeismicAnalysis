% **********************************************************************
% *                                                                    *
% *   Matlab Eigen Command                                             *
% *                                                                    *
% *   Application to Lecture 2 - Example 4                             *
% *   Version 1.0/RJA/1-10-2022                                                                 *
% *                                                                    *
% **********************************************************************

% Eigen Analysis of Lecture 2 Example 4
  ptitle =  ' Lecture 2 - Example 4'
  
% Column Properties
  h = 16*12;        % [in] 
  EI = 10000*100;   % [k-in^2]  
   
% Node Stiffnesses
 
  k11  = 60*EI/(h^3) ;
  k12  = 3*EI/(h^3) ;
  k21  = k12 ;
  k22  = 6*EI/h ;
  
  K  = [   k11,    k12; ...
           k21,     k22];  
          
% The Mass Matrix
 m  = 0.025 ;  % [k-s^2/in]
 M11 = 4*m ;
 M22 = 5*m*h^2/48 ;
 
 m = [M11,  0; ...
       0, M22];
       
fprintf('The Mass Matrix:\r'); m      

%Find the eigenpairs 
%(columns of V are eigenvectors and
% diagonals of D are eigenvalues)


fprintf('The reduced stiffness matrix:\r'); K

[V,D]=eig(K,m) ;

%Extract eigenvectors
vec1 = V(:,1) ;
vec2 = V(:,2) ;

EV = [vec1,vec2];

fprintf('The e-vectors:\r'); EV
fprintf('The e-vales:\r'); D

mu1 = D(1,1);
mu2 = D(2,2);

pi = 3.14159;
% t = (2*pi)/(D(n,n)^0.5);
t1 = (2*pi)/((mu1)^0.5);
t2 = (2*pi)/((mu2)^0.5);
w1 = (mu1)^0.5;
w2 = (mu2)^0.5;

fprintf('The periods:\r'); t1, t2
fprintf('The circular frequencies:\r'); w1, w2

% 	The Rayleigh Damping Matrix
z1 = 0.05 ; 
z2 = 0.05 ; 

freq = [ 1/w1 w1 ; 1/w2 w2 ] ;
zeta = 2 .* [ z1 ; z2 ] ;
a = inv(freq) * zeta ;
a0 = a(1);
a1 = a(2);

C = a0 .* m + a1 .* K ;

%  Point Damping
C22 = 2 * -h/4 * -h/4;
C(2,2) = C(2,2) + C22;

fprintf('The Damping Matrix:\r'); C




