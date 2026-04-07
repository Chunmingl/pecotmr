// lassosum_rss.cpp — Coordinate descent for LASSO on RSS objective
// Ported from lassosum (Mak et al 2017) functions.cpp elnet()/repelnet()
//
// Objective: min_beta  beta'R beta - 2 beta'z + 2 lambda ||beta||_1
// where R is a (possibly pre-shrunk) LD matrix and z = bhat / sqrt(n).

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Single-block coordinate descent — mirrors lassosum elnet()
// Tracks Rbeta = R * beta as a running sum (analogous to yhat = X * beta
// in the genotype version). This avoids recomputing the full dot product
// each iteration and matches the original algorithm exactly.
static int elnet_rss(double lambda1, const arma::vec& diag_R,
                     const arma::mat& R, const arma::vec& z,
                     double thr, arma::vec& beta, arma::vec& Rbeta,
                     int maxiter) {
  int p = z.n_elem;
  double dlx, del, t, bj;

  int conv = 0;
  for (int k = 0; k < maxiter; k++) {
    dlx = 0.0;
    for (int j = 0; j < p; j++) {
      bj = beta(j);
      // t = z(j) - sum_{k != j} R[j,k]*beta[k]
      //   = z(j) - (Rbeta(j) - R[j,j]*beta(j))
      t = z(j) - Rbeta(j) + diag_R(j) * bj;
      // Soft-thresholding (same as lassosum line 438)
      beta(j) = 0.0;
      if (std::abs(t) - lambda1 > 0.0)
        beta(j) = std::copysign(std::abs(t) - lambda1, t) / diag_R(j);
      if (beta(j) == bj) continue;
      del = beta(j) - bj;
      dlx = std::max(dlx, std::abs(del));
      // Update running Rbeta (analogous to yhat += del * X.col(j))
      Rbeta += del * R.col(j);
    }
    Rcpp::checkUserInterrupt();
    if (dlx < thr) {
      conv = 1;
      break;
    }
  }
  return conv;
}

// [[Rcpp::export]]
Rcpp::List lassosum_rss_rcpp(const arma::vec& z,
                              const Rcpp::List& LD,
                              const arma::vec& lambda,
                              double thr,
                              int maxiter) {
  // Compute total number of variants across all blocks
  int n_blocks = LD.size();
  std::vector<int> block_start(n_blocks), block_end(n_blocks);
  int p = 0;
  for (int b = 0; b < n_blocks; b++) {
    arma::mat Rb = Rcpp::as<arma::mat>(LD[b]);
    block_start[b] = p;
    p += Rb.n_rows;
    block_end[b] = p - 1;
  }

  if ((int)z.n_elem != p)
    Rcpp::stop("Length of z must equal total rows across all LD blocks.");

  int nlambda = lambda.n_elem;
  arma::mat beta_mat(p, nlambda, arma::fill::zeros);
  arma::ivec conv_vec(nlambda, arma::fill::zeros);
  arma::vec loss_vec(nlambda, arma::fill::zeros);
  arma::vec fbeta_vec(nlambda, arma::fill::zeros);

  // Working beta vector — warm-started across lambda path
  arma::vec beta(p, arma::fill::zeros);

  for (int i = 0; i < nlambda; i++) {
    // Block-wise coordinate descent — mirrors lassosum repelnet()
    int out = 1;
    for (int b = 0; b < n_blocks; b++) {
      arma::mat Rb = Rcpp::as<arma::mat>(LD[b]);
      int s = block_start[b];
      int e = block_end[b];
      arma::vec diag_R = Rb.diag();
      arma::vec z_blk = z.subvec(s, e);
      arma::vec beta_blk = beta.subvec(s, e);
      arma::vec Rbeta_blk = Rb * beta_blk;

      int conv_blk = elnet_rss(lambda(i), diag_R, Rb, z_blk,
                                thr, beta_blk, Rbeta_blk, maxiter);
      beta.subvec(s, e) = beta_blk;
      out = std::min(out, conv_blk);
    }

    beta_mat.col(i) = beta;
    conv_vec(i) = out;

    // Compute loss = beta'R beta - 2 z'beta (block-wise)
    double loss = -2.0 * arma::dot(z, beta);
    for (int b = 0; b < n_blocks; b++) {
      arma::mat Rb = Rcpp::as<arma::mat>(LD[b]);
      int s = block_start[b];
      int e = block_end[b];
      arma::vec beta_blk = beta.subvec(s, e);
      loss += arma::as_scalar(beta_blk.t() * Rb * beta_blk);
    }
    loss_vec(i) = loss;
    fbeta_vec(i) = loss + 2.0 * lambda(i) * arma::sum(arma::abs(beta));
  }

  return List::create(
    Named("beta") = beta_mat,
    Named("lambda") = lambda,
    Named("conv") = conv_vec,
    Named("loss") = loss_vec,
    Named("fbeta") = fbeta_vec
  );
}
