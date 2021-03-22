
// The input data is a vector 'y' of length 'N'.
data {
  real sigmaB;
  int N;
  int X;
  int J;
  matrix[X,J] f;
  matrix[X,J] d;
}
transformed data {
  // Retweet structure: Joint norm > mujx
  real E = 0;
  real D = 0;
  real F = 0;
    
  for (x in 1:X) {
    
    for (j in 1:J)
    {
      E += log(f[x, j]+1)*log(d[x,j]+1);
      E += log(d[x,j]+1);
      E += log(f[x, j]+1);
    }
  }
  //-----------------------------------------

}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  // Retweet structure: Joint norm > mujx
  real beta0;
  real betaF;
  real betad;
  real sigmaSb;
  real<lower=0> sigma;
  matrix[X,J] b;
  //-------------------------
}

transformed parameters {
  // Retweet structure: Joint norm > mujx
  int N1 = N + sigmaSb * sigmaB^(-2);
  D2 = sigmaSb * sigmaB^(-2);
  F2 = sigmaSb * sigmaB^(-2);
  Y0 = 0;
  Yf = 0;
  Yd = sigmaSb * sigmaB^(-2);
  cov_matrix[3] C;
  vector[3] mu;
  
  for (x in 1:X) {
    for (j in 1:J)
    {
      D2 += log(d[x,j]+1)^2  ;
      F2 += log(f[x,j]+1)^2 ;
      Y0 += log(b[x,j]+1) ;
      Yf += log(f[x,j]+1)*log(f[x,j]+1) ;
      Yd += log(b[x, j]+1)^2 * log(d[x,j]+1);
    }
  }
  C[1,1] = N1;
  C[2,2] = F2;
  C[3,3] = D2;
  C[1,2] = F;
  C[1,3] = D;
  C[2,3]= E;
  C = C * sigmaSb;
  C = inverse_spd(C);
  vector[3] Y;
  Y[1]= Y0;
  Y[2] = Yf;
  Y[3] = Yd;
  mu = C * Y;
  //-----------------------------------
  
  //-----------------------------------------

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
// Retweet structure: Joint norm > mujx

}

