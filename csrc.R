# ---- Sampling ----

# ---- Zi ----
# P(Zi = 1) = 0.4
Z.prob = 0.4


library(boot)
# ---- Xi|Zi ----
# P(Xi = 1|Zi = zi; gam0, gam1) = expit(gam0 + gam1 * zi) = inv.logit(gam0 + gam1 * zi) <a function in boot>

# gam0 = 0: As of Oct 1, 2018

# gam0 = logit(0.4): As of Oct 3, 2018
gam0 = 0

gam1 = log(1.25) # In R, log is natural logarithm. ( P(X = 1|Z = 0) = 2/5, P(X = 1|Z = 1) = 5/11 )

# pX.Z: value of X (0 or 1), value of Z (0 or 1), g0 (real), g1 (real) -> probability of X.Z

pX.Z = function(X, Z, g0, g1) { return( dbinom( X, size = 1, prob = inv.logit( g0 + g1 * Z ) ) ) } 


# ---- So keep an eye on X|Z ----

# ---- Wi|Xi, Zi ----
# P(Wi = 1|Xi = xi, Zi = zi) = expit(alf0 + alf1 * xi + alf2 * zi)

alf0 = logit(0.15)

alf1 = logit(0.85) - logit(0.15)

alf2 = 0

# pW.XZ: value of W (0 or 1), value of X (0 or 1), value of Z (0 or 1), 
# a0 (real), a1 (real), a2 (real) -> probability of W.XZ

pW.XZ = function(W, X, Z, a0, a1, a2) { return( dbinom( W, 
                                                        size = 1, 
                                                        prob = inv.logit( a0 + a1 * X + a2 * Z ) ) ) }


# ---- So keep an eye on W|X, Z ----

# ---- Ti|Xi, Zi ---- 
lam = 1
kap = 2
bet1 = log(0.8)
bet2 = log(1.1)

# FT.XZ: value of Time (positive), value of X (0 or 1), value of Z (0 or 1), 
# l (positive), k (positive), b1 (real), b2 (real) -> CUMULATIVE density of T.XZ

FT.XZ = function( Time, X, Z, l, k, b1, b2 ) {
  
  scl = ( l * exp(b1 * X + b2 * Z) ) ^ ( -1/k )
  
  return( pweibull(Time, shape = k, scale = scl) )
}


# FT: value of Time (positive), 
# l (positive), k (positive), b1 (real), b2 (real), g0 (real), g1 (real), pZ (between 0 and 1, P(Z = 1)) 
# -> marginal cdf of T

# We will do it by 4 partitions: (X = 0, Z = 0), (X = 0, Z = 1), (X = 1, Z = 0) and (X = 1, Z = 1)

FT = function( Time, l, k, b1, b2, g0, g1, pZ ) {
  
  FT.00 = FT.XZ(Time, 0, 0, l, k, b1, b2) * pX.Z(0, 0, g0, g1) * (1 - pZ)
  
  FT.10 = FT.XZ(Time, 1, 0, l, k, b1, b2) * pX.Z(1, 0, g0, g1) * (1 - pZ) 
  
  FT.01 = FT.XZ(Time, 0, 1, l, k, b1, b2) * pX.Z(0, 1, g0, g1) * pZ
  
  FT.11 = FT.XZ(Time, 1, 1, l, k, b1, b2) * pX.Z(1, 1, g0, g1) * pZ
  
  return( FT.00 + FT.10 + FT.01 + FT.11 )
}


# ---- Ai ----

# f.A(a; rho) = 1/rho e^(-a/rho), a > 0.

rho = 1.1

# pY.AXZ: value of Y (0 or 1), value of A (positive), value of X (0 or 1), value of Z (0 or 1)
# l (positive), k (positive), b1 (real), b2 (real) -> probability of Y.AXZ

pY.AXZ = function( Y, A, X, Z, l, k, b1, b2 ) {
  
  p.y = FT.XZ(A, X, Z, l, k, b1, b2)
  
  return( dbinom( Y, size = 1, prob = p.y ) )
}


# If one wants a more complex distribution (pdf) of A, here we go:
# pA: value of a (positive), value of mu (positive, E(A)), value of phi (positive Var(A)) -> a positive number

pA = function( A, mu, phi ) {
  
  sc = phi / mu # sc = scale
  
  sh = mu / sc # sh = shape
  
  return( dgamma( A, shape = sh, scale = sc ) )
}


# rA: value of n (natural number denoting sample size), mu, phi -> n gamma numbers

rA = function( n, mu, phi ){
  
  sc = phi / mu # sc = scale
  
  sh = mu / sc # sh = shape
  
  return( rgamma( n, shape = sh, scale = sc ) )
}

phi = rho^2 * 1.5

# mu.solve: l (positive), k (positive), b1 (real), b2 (real), g0 (real), g1 (real), pZ (aka P(Z = 1) ),
# mu.guess (positive), phi (positive), ptp (probability of testing positive of the disease), 
# tol (a minute positive number), maxit (maximum number of iteration)
# -> a real number such that P(T < A; l, k, b1, b2, g0, g1, mu, phi) = 0.5

# I need to try something first.

mu.solve = function( l, k, b1, b2, g0, g1, pZ, mu.guess, phi, ptp = 0.5, tol = 1e-6, maxit = 128 ) {
  
  # These are indicators: 0 for below (dif < -tol), and 1 (dif > tol) for over
  incre = 0.1
  
  is.over = 0
  
  # These are the guesses.
  
  guess = mu.guess
  
  obj.f = function(mu) {
    
    integrand = function(A) { return( FT(A, l, k, b1, b2, g0, g1, pZ) * pA(A, mu, phi) ) }
    
    objective = integrate( integrand, 0, Inf )$value - ptp
  }
  
  dif = obj.f( guess )
  
  # Break the loop either when we get the value within tolerance, or when we run out of iterations.
  
  niter = 0
  
  while ( abs(dif) > tol & niter < maxit) {
    
    if ( dif > tol ){ # This implies that is.over should be 1 right now.
      
      if ( is.over == 0 ) { incre = incre / 2 } 
      
      # However, it is 0 before; a change in sign requires decreasing increment.
      
      is.over = 1 # Because dif > tol, is.over needs resetting to 1, regardless.
      
      guess = guess - incre # Because dif > tol, I need to decrease the guess, regardless. The difference is how much.
      
      dif = obj.f( guess )
      
    } else if ( dif < -tol ) {
      
      if ( is.over == 1 ) { incre = incre / 2} 
      
      is.over = 0 # Because dif < -tol, is.over needs resetting to 0, regardless.
      
      guess = guess + incre # Because dif > tol, I need to decrease the guess, regardless.
      
      dif = obj.f( guess )
      
    }
    
    niter = niter + 1
  }
  
  if ( abs(dif) <= tol ) { return(guess) } else {
    
    warning("You have reached the maximum number of iterations, but the algorithm did not converge.")
    
    return(guess)
  }
  
}

mu = mu.solve( lam, kap, bet1, bet2, gam0, gam1, Z.prob, mu.guess = 1, phi, ptp = 0.5)

# E(A) = 1.258, Var(A) = 1.815

# ---- Ri ----
# P(Ri = 1) = 0.3 so we have 70% missing data. 
p.r = 0.3

# ---- Preparation ----

# pX.YAWZ: value of X (0 or 1), value of Y (0 or 1), value of A (positive), value of W (0 or 1), value of Z (0 or 1),
# l (positive), k (positive), b1 (real), b2 (real), a0 (real), a1 (real), a2 (real), g0 (real), g1 (real)

pX.YAWZ = function( X, Y, A, W, Z, l, k, b1, b2, a0, a1, a2, g0, g1 ) {
  
  part = function(X) { 
    
    pY.AXZ(Y, A, X, Z, l, k, b1, b2) * pW.XZ(W, X, Z, a0, a1, a2) * pX.Z(X, Z, g0, g1) }
  
  return( part(X) / ( part(0) + part(1) ) )
} 


# poids: (entire) dataframe of current status (R, Y, A, X, W, Z), l (positive), k (positive), 
# b1 (real), b2 (real), a0 (real), a1 (real), a2 (real), g0 (real), g1 (real) -> a vector of weights
# If R = 1, the weight must be 1; otherwise, it will be P(X = x|Y = y, A = a, W = w, Z = z).

poids = function(cs, l, k, b1, b2, a0, a1, a2, g0, g1) {
  
  les.poids = as.numeric()
  
  for ( h in 1 : nrow( cs ) ) {
    
    dline.h = cs[ h, ]
    
    if (dline.h$R == 0) {
      
      X.h = dline.h$X
      
      Y.h = dline.h$Y
      
      A.h = dline.h$A
      
      W.h = dline.h$W
      
      Z.h = dline.h$Z
      
      les.poids[h] = pX.YAWZ(X.h, Y.h, A.h, W.h, Z.h, l, k, b1, b2, a0, a1, a2, g0, g1)
    } else { les.poids[h] = 1 }
  }
  
  return( les.poids )
  
}


# ---- Score functions (Note functional independece) ----

# score.Y: Y (0 or 1), A (positive), X (0 or 1), Z (0 or 1), 
# l (positive), k (positive), b1 (real), b2 (real)

score.Y = function( Y, A, X, Z, l, k, b1, b2 ) {
  
  bets = c( log(l), b1, b2, k )
  
  covar = c( 1, X, Z, log(A) )
  
  mult = exp( sum( covar * bets ) )
  
  p = 1 - exp( -mult )
  
  return ( as.vector( covar * mult * ( Y / p - 1 ) ) )
}


# score.W: W (0 or 1), X (0 or 1), Z (0 or 1), a0 (real), a1 (real), a2 (real) 
# -> a real vector in R^3

score.W = function( W, X, Z, a0, a1, a2 ) {
  
  mult = W - inv.logit( a0 + a1 * X + a2 * Z )
  
  covar = c( 1, X, Z )
  
  return( covar * mult )
}


# score.X: X (0 or 1), Z (0 or 1), g0 (real), g1 (real) -> a real vector in R^2

score.X = function( X, Z, g0, g1 ) {
  
  mult = X - inv.logit( g0 + g1 * Z )
  
  covar = c( 1, Z )
  
  return( covar * mult )
}


# sc.sct: binary (X), binary (Y), positive-continuous (A), binary (W), binary (Z),
# positive-continuous (l), positive-continuous (k), continuous (b1), continuous (b2),
# continuous (a0), continuous (a1), continuous (a2), continuous (g0), continuous (g1)
# -> a 9*9 symmetric matrix: score times score transposed

sc.sct = function( X, Y, A, W, Z, l, k, b1, b2, a0, a1, a2, g0, g1 ) {
  
  sc.overall = c( score.Y(Y, A, X, Z, l, k, b1, b2), score.W(W, X, Z, a0, a1, a2),
                  score.X(X, Z, g0, g1) )
  
  return( sc.overall %*% t( sc.overall ) )
}


# ---- About the data ----

# augmenter: a dataframe of current status -> an augmented dataframe of current status
# (When R = 1, the dataline stays as is; when R = 0, there will be 2 replicates, one for
#  when X = 1 and the other for when X = 0.)

augmenter = function( cs ) {
  
  cs.complete = cs[ cs$R == 1,]
  
  if ( nrow( cs.complete ) == nrow(cs) ) {
    
    # If the complete dataset has the same number of line (observation) as the original set,
    # That means there's nothing missing.
    
    return( cs ) } else {
      
      row.names( cs.complete ) = 1 : nrow( cs.complete )
      
      cs.missing = cs[ cs$R == 0, ]
      
      row.names( cs.missing ) = ( nrow( cs.complete ) + 1 ) : nrow( cs )
      
      # Blind fold:
      cs.missing$X = NA
      
      # Copy:
      csm.copy = cs.missing
      
      row.names( csm.copy ) = ( nrow(cs) + 1 ) : ( nrow(cs) + nrow(csm.copy) )
      
      # Assign X = 0 to cs.missing and X = 1 to csm.copy:
      cs.missing$X = 0
      
      csm.copy$X = 1
      
      cs.augmented = rbind( cs.complete, cs.missing, csm.copy )
      
      return( cs.augmented )
    }
  
}


# esst: a dataframe of current status (R, Y, A, X, W, Z), l (positive), k (positive), 
# b1 (real), b2 (real), a0 (real), a1 (real), a2 (real), g0 (real), g1 (real)
# -> a 9 * 9 symmetric real matrix (conditional expectation of score * score transposed)
# (We need to use the augmented version of dataframe.)

esst = function( cs.augmented, l, k, b1, b2, a0, a1, a2, g0, g1 ) {
  
  cs.incomplete = cs.augmented[ cs.augmented$R == 0, ]
  
  poids.manquants = poids( cs.incomplete, l, k, b1, b2, a0, a1, a2, g0, g1 )
  
  condvar = matrix( rep(0, 81), nrow = 9, byrow = T )
  
  n.incomplete = nrow( cs.incomplete ) / 2 # How many missing data are actually out there? Half.
  
  for ( h in 1 : n.incomplete ) {
    
    dline.h0 = cs.incomplete[ h, ]
    
    X.h0 = dline.h0$X
    
    Y.h0 = dline.h0$Y
    
    A.h0 = dline.h0$A
    
    W.h0 = dline.h0$W
    
    Z.h0 = dline.h0$Z
    
    poids.h0 = poids.manquants[h]
    
    dline.h1 = cs.incomplete[ h + n.incomplete, ]
    
    X.h1 = dline.h1$X
    
    Y.h1 = dline.h1$Y
    
    A.h1 = dline.h1$A
    
    W.h1 = dline.h1$W
    
    Z.h1 = dline.h1$Z
    
    poids.h1 = poids.manquants[ h + n.incomplete ]
    
    # 2nd-moment contribution from h
    sm.h = sc.sct( X.h0, Y.h0, A.h0, W.h0, Z.h0, l, k, b1, b2, a0, a1, a2, g0, g1 ) * poids.h0 +
           sc.sct( X.h1, Y.h1, A.h1, W.h1, Z.h1, l, k, b1, b2, a0, a1, a2, g0, g1 ) * poids.h1
    
    # squared-1st-moment contribution from h
    fm.h = c( score.Y( Y.h0, A.h0, X.h0, Z.h0,l, k, b1, b2 ), 
              score.W( W.h0, X.h0, Z.h0, a0, a1, a2), 
              score.X( X.h0, Z.h0, g0, g1 ) ) * poids.h0 + 
      
           c( score.Y( Y.h1, A.h1, X.h1, Z.h1,l, k, b1, b2 ), 
              score.W( W.h1, X.h1, Z.h1, a0, a1, a2), 
              score.X( X.h1, Z.h1, g0, g1 ) ) * poids.h1
    
    cvar.h = sm.h - fm.h %*% t( fm.h )
    
    condvar = condvar + cvar.h
  }
  
  return( condvar )
}

# ---- EM Algorithm. ---- 

# em.cs: dataframe of current status (not yet augmented), 
# niter (number of iterations, default = 12800), tol (tolerance, default = 1e-6) 
# -> listof estimates of bet0 (ln lam), bet1, bet2, bet3 (kap), alf0, alf1, alf2, gam0, gam1,
#       and EM covariance matrix

em.cs = function(cs, niter = 12800, tol = 1e-6) {
  
  cs.pseudo = augmenter( cs )
  
  # Initialize estimates based on 'augmented' data (initial weights = 1)
  Ysum.init = summary( glm( Y ~ X + Z + log(A), family = binomial(link = 'cloglog'), data = cs.pseudo ) )
  
  bets.init = Ysum.init$coefficients[1:4]
  
  Wsum.init = summary( glm( W ~ X + Z, family = binomial(link = 'logit'), data = cs.pseudo ) )
  
  alfs.init = Wsum.init$coefficients[1:3]
  
  Xsum.init = summary( glm( X ~ Z, family = binomial(link = 'logit'), data = cs.pseudo ) )
  
  gams.init = Xsum.init$coefficients[1:2]
  
  current.estimate = c( bets.init, alfs.init, gams.init ) 
  
  for (h in 1:niter) {
    
    poids.h = poids(cs.pseudo, l = exp(current.estimate[1]), k = current.estimate[4], 
                    b1 = current.estimate[2], b2 = current.estimate[3], 
                    a0 = current.estimate[5], a1 = current.estimate[6], 
                    a2 = current.estimate[7], g0 = current.estimate[8], 
                    g1 = current.estimate[9])
    
    # Note: We use weights at step h to get estimates at step h+1.
    
    Ysum.h = summary( glm( Y ~ X + Z + log(A), family = binomial(link = 'cloglog'), 
                           data = cs.pseudo, weights = poids.h ) )
    
    bets.h = Ysum.h$coefficients[1:4]
    
    Wsum.h = summary( glm( W ~ X + Z, family = binomial(link = 'logit'), 
                           data = cs.pseudo, weights = poids.h ) )
    
    alfs.h = Wsum.h$coefficients[1:3]
    
    Xsum.h = summary( glm( X ~ Z, family = binomial(link = 'logit'), 
                           data = cs.pseudo, weights = poids.h ) )
    
    gams.h = Xsum.h$coefficients[1:2]
    
    estimate.h = c( bets.h, alfs.h, gams.h )
    
    dif.h = norm( estimate.h - current.estimate, type = '2')
    
    if (dif.h < tol) { # tolerance satisfied
      
      # First expression in calculating EM Info (B.matrix)
      YHess.final = Ysum.h$cov.unscaled
      
      WHess.final = Wsum.h$cov.unscaled
      
      XHess.final = Xsum.h$cov.unscaled
      
      covar = rbind( cbind( YHess.final, matrix( rep(0, 20), nrow = 4, byrow = T ) ),
                        
                     cbind( matrix( rep(0, 12), nrow = 3, byrow = T ), WHess.final,
                            matrix( rep(0, 6), nrow = 3, byrow = T ) ),
                        
                     cbind( matrix( rep(0, 14), nrow = 2, byrow = T ), XHess.final ) )
      
      covar = unname( covar )
      
      B.matrix = solve( covar )
      
      # Second expression in calculating EM Info ( conditional weighted average of score * t(score) )
      # Recalculating weights
      l.final = exp( bets.h[1] )
      
      k.final = bets.h[4]
      
      b1.final = bets.h[2]
      
      b2.final = bets.h[3]
      
      a0.final = alfs.h[1]
      
      a1.final = alfs.h[2]
      
      a2.final = alfs.h[3]
      
      g0.final = gams.h[1]
      
      g1.final = gams.h[2]

      S2.matrix = esst( cs.pseudo, l.final, k.final, b1.final, b2.final, a0.final, a1.final, a2.final,
                        g0.final, g1.final )
      
      Info = B.matrix - S2.matrix
      
      EM.VAR = solve( Info )
      
      output = list( estimate.h, EM.VAR )
      
      return( output )
      
    } else if ( h < niter ) { 
        
      current.estimate = estimate.h 
      
    } else {
      
      # First expression in calculating EM Variance
      YHess.final = Ysum.h$cov.unscaled
      
      WHess.final = Wsum.h$cov.unscaled
      
      XHess.final = Xsum.h$cov.unscaled
      
      covar = rbind( cbind( YHess.final, matrix( rep(0, 20), nrow = 4, byrow = T ) ),
                        
                     cbind( matrix( rep(0, 12), nrow = 3, byrow = T ), WHess.final,
                            matrix( rep(0, 6), nrow = 3, byrow = T ) ),
                        
                     cbind( matrix( rep(0, 14), nrow = 2, byrow = T ), XHess.final ) )
      
      covar = unname( covar )
      
      B.matrix = solve( covar )
      
      # Recalculating weights
      l.final = exp( bets.h[1] )
      
      k.final = bets.h[4]
      
      b1.final = bets.h[2]
      
      b2.final = bets.h[3]
      
      a0.final = alfs.h[1]
      
      a1.final = alfs.h[2]
      
      a2.final = alfs.h[3]
      
      g0.final = gams.h[1]
      
      g1.final = gams.h[2]
      
      S2.matrix = esst( cs.pseudo, l.final, k.final, b1.final, b2.final, a0.final, a1.final, a2.final,
                        g0.final, g1.final )
      
      Info = B.matrix - S2.matrix
      
      EM.VAR = solve( Info )
      
      output = list( estimate.h, EM.VAR )
      
      print('You have reached the maximum number of iterations, but the algorithm did not converge.')
      
      return( output )
      
    }
    
  }
}

# true values put in the vector form

tru.val = c( log(lam), bet1, bet2, kap, alf0, alf1, alf2, gam0, gam1 )

# ---- Now, reproduce 2000 times ----  

# reproduce: nset (number of dataset, default 2000), nline (number of observation per set, default 500)
# pZ P(Z = 1) (default 0.4), b1 (real, default bet1), b2 (real, default bet1), 
# a0 (real, default alf0), a1 (real, default alf1), a2 (real, default alf2), 
# g0 (real, default gam0), g1 (real, default gam1), l (positive, default lam), k (positive, default kappa), 
# mu (positive, default mu), phi (positive, default phi), pR P(R = 1) (default 0.3)

# -> a 27-element vector that contains estimates, standard errors, and coverage probabilities
# (9 parameter, 3 pieces of summary - estimate, standard error, and coverage probability each)

reproduce = function( nline = 500, pZ = Z.prob, b1 = bet1, b2 = bet2, 
                      
                      a0 = alf0, a1 = alf1, a2 = alf2, 
                      
                      g0 = gam0, g1 = gam1, l = lam, k = kap, p = rho, pR = p.r) {
  
  # Z:
  Z.vec = rbinom( n = nline, size = 1, prob = pZ )
  
  # X|Z
  X.vec = as.numeric()
  
  for ( i in 1 : nline ) {
    
    p.i = pX.Z( X = 1, Z = Z.vec[i], g0, g1 )
    
    X = rbinom( n = 1, size = 1, prob = p.i )
    
    X.vec[i] = X
  }
  
  # W|X, Z
  W.vec = as.numeric()
  
  for ( j in 1 : nline ) {
    
    p.j = pW.XZ( W = 1, X = X.vec[j], Z = Z.vec[j], a0, a1, a2 )
    
    W = rbinom( n = 1, size = 1, prob = p.j )
    
    W.vec[j] = W
  }
  
  # T|X, Z
  T.vec = as.numeric()
  
  for ( m in 1 : nline ) {
    
    scl.m = ( l * exp( b1 * X.vec[m] + b2 * Z.vec[m] ) ) ^ ( -1/k )
    
    Time = rweibull( n = 1, shape = k, scale = scl.m )
    
    T.vec[m] = Time
  }
  
  # A:
  A.vec = rA( n = nline, mu = mu, phi = phi )
  
  # Y:
  Y.vec = as.numeric( T.vec <= A.vec )
  
  # R:
  R.vec = rbinom( n = nline, size = 1, prob = pR )
  
  # Forming dataframe:
  current.status = data.frame( cbind( R.vec, Y.vec, A.vec, X.vec, W.vec, Z.vec ) )
  
  colnames( current.status ) = c( 'R', 'Y', 'A', 'X', 'W', 'Z' )
  
  # Now EM:
  em.output = em.cs( current.status )
  
  estimates = em.output[[1]]
  
  em.variance = em.output[[2]] # This is a matrix
  
  # Summary:
  esc = 1:27
  
  for ( o in 1:9 ) {
    
    est.o = estimates[o]
    
    se.o = sqrt( em.variance[o,o] )
    
    tru.o = tru.val[o]
    
    lb.o = est.o - 1.96 * se.o
    
    ub.o = est.o + 1.96 * se.o
    
    cover.o = as.numeric( ( lb.o <= tru.o ) & ( tru.o <= ub.o ) )
    
    esc[ 3 * o - 2 ] = est.o
    
    esc[ 3 * o - 1 ] = se.o
    
    esc[ 3 * o ] = cover.o
  }
  
  return( esc )
}


# em.summary: est.xx (2000*18 matrix) -> a 9*5 matrix, detailing average of estimates, standard error, 
# true value, bias, and coverage probability of 95% confidence interval for each parameter.

em.summary = function(est.xx) {
  
  tab = matrix( NA, nrow = 9, ncol = 6, byrow = T )
  
  colnames( tab ) = c( 'average of estimates', 'sample standard error', 'Louis standard error', 'true value', 'bias',
                        'coverage probability of nominal 95% confidence interval' )
  
  row.names( tab ) = c( 'ln(lamda)', 'beta1', 'beta2', 'kappa', 'alpha0', 'alpha1', 'alpha2',
                         'gamma0', 'gamma1' )
  
  for (i in 1:9) {
    
    estimate.i = mean( est.xx[, 3 * i - 2 ] )
    
    smpse.i = sd( est.xx[, 3 * i - 2 ] )
    
    louis.i = mean( est.xx[, 3 * i - 1 ] )
    
    tru.i = tru.val[i]
    
    bias.i = estimate.i - tru.i
    
    cover.i = mean( est.xx[, 3 * i] )
    
    tab[i,] = c(estimate.i, smpse.i, louis.i, tru.i, bias.i, cover.i)
    
  }
  
  return( tab )
}

