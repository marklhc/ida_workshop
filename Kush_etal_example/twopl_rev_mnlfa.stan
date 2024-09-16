//
// From https://mc-stan.org/docs/stan-users-guide/regression.html#item-response-models.section
//

// The input data is an array of integers.
data {
  int<lower=1> J;                     // number of students
  int<lower=1> K;                     // number of questions
  int<lower=1> N;                     // number of observations
  array[N] int<lower=1, upper=J> jj;  // student for observation n
  array[N] int<lower=1, upper=K> kk;  // question for observation n
  array[N] int<lower=0, upper=1> y;   // correctness for observation n
  int<lower=1> p;                     // number of covariates
  matrix[N, p] X;                     // covariate predicting latent variance
}
parameters {
  vector[J] alpha_star;      // ability of student j - mean ability
  vector[K] gamma0_star;     // discrimination of question k (before rotation)
  vector[K] beta0;           // difficulty of question k
  vector[p] g_star;          // effect of x on latent mean
  vector[p] h;               // effect of x on log(latent SD)
  // Horseshoe prior for paths on intercepts
  matrix[p, K] zb;           // unregularized effect of x on intercepts
  real<lower=0> taub;        // global shrinkage
  matrix<lower=0>[p, K] lb;  // local shrinkage
  // Horseshoe prior for paths on loadings
  matrix[p, K] zc_star;      // unregularized effect of x on loadings/discrimination
  real<lower=0> tauc;        // global shrinkage
  matrix<lower=0>[p, K] lc;  // local shrinkage
}
transformed parameters {
  matrix[p, K] b = zb .* lb * taub;
  matrix[p, K] c_star = zc_star .* lc * tauc;
}
// The model to be estimated. 
model {
//   alpha_star ~ std_normal();         // informative true prior
  h ~ std_normal();
  g_star ~ std_normal();
  alpha_star ~ std_normal();
  beta0 ~ normal(0, 2.5);          // informative true prior
  gamma0_star ~ normal(1, 2.5);      // informative true prior
  to_vector(zb) ~ std_normal();
  to_vector(lb) ~ student_t(3, 0, 1);
  taub ~ std_normal();
  to_vector(zc_star) ~ std_normal();
  to_vector(lc) ~ student_t(3, 0, 1);
  tauc ~ std_normal();
  {
    vector[N] dbeta;  // deviation from reference intercepts
    vector[N] dgamma;  // deviation from reference intercepts
    for (n in 1:N) {
      dbeta[n] = X[n, ] * b[, kk[n]];  // choose the right item for observation n;
      dgamma[n] = X[n, ] * c_star[, kk[n]];
    }
    y ~ bernoulli_logit((gamma0_star[kk] + dgamma) .*
                            (exp(X * h) .* (alpha_star[jj] + X * g_star)) +
                            beta0[kk] + dbeta);
  }
}
generated quantities {
    vector[K] gamma0;
    vector[J] alpha;
    vector[p] g;
    matrix[p, K] c;
    {
        int sign_l1 = gamma0_star[1] > 0 ? 1 : -1;
        gamma0 = sign_l1 * gamma0_star;
        alpha = sign_l1 * alpha_star;
        g = sign_l1 * g_star;
        c = sign_l1 * c_star;
    }
}
