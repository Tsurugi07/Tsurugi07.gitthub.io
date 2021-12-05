A1 = [3 -0.1 -0.2; 0.1 7 -0.3; 0.3 -0.2 10] # Création d'une matrice A1 

B1 = [7.85; -19.3;71.4] # Création d'un vecteur B1 

function Jacobi_ite(A1,B1)
  (n,p)= size(A1); # Dénifir la taille de la Matrice A1
  X1 = zeros(n,1); # Affecter la taille pour la matrice X1
  for k = 1:11 # Jusqu'a 11 itérations
    for i = 1 : n 
      # Composante du vecteur X
     X[i] = (B1[i] - A1[i,[1:i-1,i+1:n]]* X1([1:i-1,i+1:n]))/A1[i,i];  
     X' # voir les iterations effectuées et la matrice finale
     k # nombre d'itérations effectuées 
    end 
    err_abso = ((X[i] - X1[i])/X[i])*100; # Calcul de l'erreur
    X1 = X';
    if err_abso <= 0.0001
      break
    end
  end
return X;
end

Jacobi_ite(A1,B1)

A3 = [3 -2 0;-1 3 -2;0 -1 3] # matrice A

B3 = [1;1;1] # vecteur B3

Xt = [17/15;6/5;11/15] # vecteur X théorique

k = 11 # nombre d'itérations

Xo = [0,0,0] # vecteur X initial

function Jacobideuxieme(A3,B3,Xo,k)
  (n,p) = size(A3); 
  Xo = zeros(n,1);
  for m = 1:k # 11 iterations
    for i = 1 : n
      Y[i] = (B3[i] - A3[i,[1:i-1,i+1:n]]* Xo([1:i-1,i+1:n]))/A3[i,i];
      Y' # Afficher chaque vecteur Y après calcul 
      err_abso = ((Y[i] - Xo[i])/Y[i])*100;# calcul de l'erreur
      err_abso # Afficher chaque erreur
      m  # décompte du nombre d'itérations 
    end
    Xo = Y';
    if err_abso <= 0.0001
      break
    end
  end
return Y;
end

Jacobideuxieme(A3,B3,Xo,k)

tol = 10^-12;# Nous devons atteindre 10^-12

A = [3 -2 0;-1 3 -2;0 -1 3];# matrice A

B = [1;1;1]; # vecteur B

X= [0;0;0];

# Debut de la fonction avec X= [0;0;0]
# Elle doit renvoyer la matrice finale ainsi que le nombre d'iérations
function jacobi(A,B,X,tol)
n = length(B);
for j = 1 : n
     # Composantes du vecteur X
     X[j] = ((B[j] - A[j,[1:j-1,j+1:n]] * X3([1:j-1,j+1:n])) / A[j,j]); 
end
x1 = X;
iter = 0; # nombre d'ireration que l'on fixe a zero 
while norm(x1-X,1) > tol # 
    for j = 1 : n
      # Calcul du vecteur n_ny
     x_ny[j] = ((B[j] - A[j,[1:j-1,j+1:n]] * x1([1:j-1,j+1:n])) / A[j,j]);
    end
    X = x1;
    x1 = x_ny'; # x1 = à la transposée 
    x = x1;
end
    iter = j; # determiner le nombre d'iterations final
    iter # Afficher le nombre d'iterations final
    return x1;
    return iter;
end

jacobi(A,B,X,tol)

# Debut de la fonction avec X2= [1;1;1]
# Elle doit renvoyer la matrice finale ainsi que le nombre d'iérations
function jacobi(A,B,X2,tol)
n = length(B);
for j = 1 : n
   # Composantes du vecteur X
     X[j] = ((B[j] - A[j,[1:j-1,j+1:n]] * X3([1:j-1,j+1:n])) / A[j,j]); 
end
x1 = X;
iter = 0; # nombre d'ireration que l'on fixe a zero 
while norm(x1-X,1) > tol # 
    for j = 1 : n
      # Calcul du vecteur n_ny
     x_ny[j] = ((B[j] - A[j,[1:j-1,j+1:n]] * x1([1:j-1,j+1:n])) / A[j,j]);
    end
    X = x1;
    x1 = x_ny'; # x1 = à la transposée 
    x = x1;
end
    iter = j; # determiner le nombre d'iterations final
    iter # Afficher le nombre d'iterations final
    return x1;
    return iter;
end

jacobi(A,B,X2,tol)

# Debut de la fonction avec X3= [1;10;100]
# Elle doit renvoyer la matrice finale ainsi que le nombre d'iérations
function jacobi(A,B,X3,tol)
n = length(B);
for j = 1 : n
   # Composantes du vecteur X
     X[j] = ((B[j] - A[j,[1:j-1,j+1:n]] * X3([1:j-1,j+1:n])) / A[j,j]); 
end
x1 = X;
iter = 0; # nombre d'ireration que l'on fixe a zero 
while norm(x1-X,1) > tol # 
    for j = 1 : n
      # Calcul du vecteur n_ny
     x_ny[j] = ((B[j] - A[j,[1:j-1,j+1:n]] * x1([1:j-1,j+1:n])) / A[j,j]);
    end
    X = x1;
    x1 = x_ny'; # x1 = à la transposée 
    x = x1';
end
    iter = j; # determiner le nombre d'iterations final
    iter # Afficher le nombre d'iterations final
    return x1;
    return iter;
end

jacobi(A,B,X3,tol)

A4 = [1 1 2;-1 0 1;2 1 2]# matrice A

B4 = [4;3;3] # vecteur B3

Xt = [17/15;6/5;11/15] # vecteur X théorique

k = 11 # nombre d'itération 

Xo = [0;0;0] # vecteur X initial 

Xo1 = [1;1;1] # secong vecteur X initial 

function Jacobitroisieme(A4,B4,Xo,k)
  (n,p) = size(A4); 
  Xo = zeros(n,1);
  for m = 1:k # 11 iterations
    for i = 1 : n
      X[i] = (B4[i] - A4[i,[1:i-1,i+1:n]]* Xo([1:i-1,i+1:n]))/A4[i,i];
      X'
    endfor
    err_abso = ((X[i] - Xo[i])/X[i])*100;
    Xo = X';
    if err_abso <= 0.0001
      break
    end 
  end
return X;
return m;
end

Jacobitroisieme(A4,B4,Xo,k)

# Effectuer le meme programme mais avec pour vecteur initial X = [1;1;1)
Xo1 = [1;1;1]; 
function Jacobitroisieme(A4,B4,Xo1,k)
  (n,p) = size(A4); #
  Xo1 = zeros(n,1);
  for m = 1:k # 11 iterations
    for i = 1 : n
      X[i] = (B4[i] - A4[i,[1:i-1,i+1:n]]* Xo1([1:i-1,i+1:n]))/A4[i,i];
      X'
    end
    err_abso = ((X[i] - Xo1[i])/X[i])*100;
    Xo1 = X';
    if err_abso <= 0.0001
      break
    end 
  end
return X;
return m;
end

Jacobitroisieme(A4,B4,Xo1,k)

n = 10;

# utiliser la fonction TRID pour obetnir matrice 10*10
function TriD(n)
  A=zeros(n,n);
  A[n,n]=3;
  for i=1:n-1
    A[i,i]=3;
    A[i+1,i]=-2;
    A[i,i+1]=-1;
  end
  return A;
end

TriD(10)

# Concevoir la fonction 
#### (Attention celle ci n'est pas complète)####
function J(A) # matrice A issue de la fonction TriD
    # A = TriD(n);
    (n,p)=size(A); 
    X = ones(n,1); 
    B = zeros(n,1);
    B[1] = X[1]*A[1,1]; # calcul de la premiere valeur de B
    for m = 1:11  
      for j = 1 : n  
        # calcul de la suite des composants du vecteur B  
           B[j] = A[j,(1:j-1,j+1:n)] * X[1:j-1,j+1:n] / A[j,j] + X[j]; 
           B # 
           m # nombre d'itérations
      end
    end
return B;
end 

J(A)

  D = diag(diag(A)); # matrice diagonale de A
  E = -triu(A)+D; # matrice superieure de A
  F = -tril(A)+D; # matrice inferieure de A


