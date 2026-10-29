// Numerical reductions for the R likelihood engine. Model equations, analytic
// derivatives, quadrature rules and gradient assembly stay in R; these helpers
// only reduce arrays that R has laid out.
//
// Conventions. Matrices are R's column-major arrays. `ptr` holds the zero-based
// cumulative row offsets of consecutive panels (length panels + 1, R's
// group_ptr). Log terms of n mixtures over m nodes arrive as one vector with the
// mixture index varying fastest (an n x m matrix). A mixture row containing NaN
// or +Inf gives a NaN value and NaN weights, so R rejects the evaluation; a row
// of -Inf (no mass) gives value -Inf and zero weights. Exactly zero weights mask
// their nodes' derivatives, so an infinite partial at an underflowed node cannot
// produce NaN; plain dot products elsewhere keep IEEE 0 * Inf = NaN.
#include <RcppArmadillo.h>
#include <climits>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

arma::uvec panel_boundaries(const arma::vec& ptr, arma::uword rows) {
    if (ptr.is_empty() || !ptr.is_finite() || arma::any(ptr != arma::floor(ptr)) ||
        ptr[0] != 0 || ptr[ptr.n_elem-1] != rows || arma::any(arma::diff(ptr) < 0))
        Rcpp::stop("Panel boundaries must be ordered integers spanning all rows");
    return arma::conv_to<arma::uvec>::from(ptr);
}

// Normalized weights (n x m) of n mixtures; their log sums are returned in value.
arma::mat normalized_weights(const arma::vec& terms, double rows, arma::vec& value) {
    if (!std::isfinite(rows) || rows < 0 || rows != std::floor(rows) || rows > INT_MAX)
        Rcpp::stop("Mixture rows must be a nonnegative integer within matrix limits");
    const arma::uword n = rows;
    if (n ? terms.n_elem % n != 0 : terms.n_elem != 0)
        Rcpp::stop("Log terms and mixture rows must agree");
    arma::mat weights = arma::reshape(terms, n, n ? terms.n_elem/n : 0);
    value.set_size(n);
    value.fill(R_NegInf);
    for (arma::uword i = 0; i < n && weights.n_cols; ++i) {
        auto row = weights.row(i);
        if (row.has_nan() || arma::any(row == R_PosInf)) {
            value[i] = R_NaN;
            row.fill(R_NaN);
        } else if (row.max() == R_NegInf) row.zeros();
        else {
            const double high = row.max();
            row = arma::exp(row-high);
            const double mass = arma::accu(row);
            row /= mass;
            value[i] = high + std::log(mass);
        }
    }
    return weights;
}

// Exactly zero mixture weights mask absent-node derivatives before multiplication.
arma::mat weighted_rows(const arma::mat& x, const arma::mat& weights) {
    if (x.n_rows != weights.n_elem) Rcpp::stop("Log terms and score rows must agree");
    arma::mat out(weights.n_rows, x.n_cols, arma::fill::zeros);
    if (x.is_empty()) return out;
    const arma::uvec absent = arma::find(weights == 0);
    for (arma::uword j = 0; j < x.n_cols; ++j) {
        arma::mat values = arma::reshape(x.col(j), weights.n_rows, weights.n_cols);
        values.elem(absent).zeros();
        out.col(j) = arma::sum(values % weights, 1);
    }
    return out;
}

// Column sums of x over each panel's rows: panels x ncol(x).
// [[Rcpp::export(rng = false)]]
arma::mat group_sum_cpp(const arma::mat& x, const arma::vec& ptr) {
    const arma::uvec ends = panel_boundaries(ptr, x.n_rows);
    arma::mat out(ends.n_elem-1, x.n_cols, arma::fill::zeros);
    if (x.is_empty()) return out;
    for (arma::uword g = 0; g < out.n_rows; ++g)
        if (ends[g+1] > ends[g]) out.row(g) = arma::sum(x.rows(ends[g], ends[g+1]-1), 0);
    return out;
}

// Panel-fast rows; the final column sums weights. Optional node weights collapse first.
// [[Rcpp::export(rng = false)]]
arma::mat group_weighted_design_cpp(const arma::mat& x, const arma::mat& weights,
        const arma::vec& ptr, const arma::mat& mixture, bool collapse) {
    const arma::uvec ends = panel_boundaries(ptr, x.n_rows);
    const arma::uword groups = ends.n_elem-1, blocks = collapse ? 1 : weights.n_cols;
    if (weights.n_rows != x.n_rows) Rcpp::stop("Design and weights must have matching rows");
    if (collapse && (mixture.n_rows != groups || mixture.n_cols != weights.n_cols))
        Rcpp::stop("Mixture weight dimensions do not agree");
    arma::mat out(groups*blocks, x.n_cols+1, arma::fill::zeros);
    // A single weight column needs no temporary panel matrices or scatter map.
    // Ordinary dot products retain 0 * Inf = NaN; mixture masking stays below.
    if (!collapse && blocks == 1) {
        for (arma::uword g = 0; g < groups; ++g) if (ends[g+1] > ends[g]) {
            const arma::uword first = ends[g], last = ends[g+1]-1;
            const auto w = weights.col(0).subvec(first, last);
            for (arma::uword j = 0; j < x.n_cols; ++j)
                out(g,j) = arma::dot(x.col(j).subvec(first, last), w);
            out(g,x.n_cols) = arma::sum(w);
        }
        return out;
    }
    for (arma::uword g = 0; g < groups && blocks; ++g) if (ends[g+1] > ends[g]) {
        arma::mat w = weights.n_cols ? arma::mat(weights.rows(ends[g], ends[g+1]-1)) :
            arma::zeros<arma::mat>(ends[g+1]-ends[g], 1);
        if (collapse && weights.n_cols) {
            const arma::vec node = mixture.row(g).t();
            const arma::uvec keep = arma::find(node != 0);
            w = keep.is_empty() ? arma::zeros<arma::mat>(w.n_rows, 1) :
                arma::mat(w.cols(keep)*node.elem(keep));
        }
        arma::mat reduced = arma::sum(w, 0).t();
        if (x.n_cols) reduced = arma::join_rows(w.t()*x.rows(ends[g], ends[g+1]-1), reduced);
        out.rows(g+groups*arma::regspace<arma::uvec>(0, blocks-1)) = reduced;
    }
    return out;
}

// Log sums of n mixtures and their normalized weights.
// [[Rcpp::export(rng = false)]]
Rcpp::List logsum_weights_cpp(const arma::vec& terms, double rows) {
    arma::vec value;
    const arma::mat weights = normalized_weights(terms, rows, value);
    return Rcpp::List::create(Rcpp::Named("value")=Rcpp::NumericVector(value.begin(), value.end()),
        Rcpp::Named("weights")=weights);
}

// Log sums of n mixtures in column 1, then mixture-weighted score columns; the
// rows of score follow the log terms.
// [[Rcpp::export(rng = false)]]
arma::mat logsum_score_cpp(const arma::vec& log_term, const arma::mat& score, double rows_input) {
    arma::vec value;
    const arma::mat weights = normalized_weights(log_term, rows_input, value);
    return arma::join_rows(value, weighted_rows(score, weights));
}

namespace {
inline void compensated_add(double value, double& sum, double& correction) {
    const double adjusted = value-correction;
    const double next = sum+adjusted;
    correction = (next-sum)-adjusted;
    sum = next;
}
}

// Generic log mass a[i]+b[j]+cross*x[i]*y[j], with axis marginals and first
// cross moments. R supplies the rule and assembles model scores. Two streaming
// passes use O(length(x)+length(y)) memory, never a tensor-sized buffer.
// [[Rcpp::export(rng = false)]]
Rcpp::List bilinear_reduce_cpp(const arma::vec& a, const arma::vec& b,
        const arma::vec& x, const arma::vec& y, double cross, bool moments = true) {
    const arma::uword nu = a.n_elem, nv = b.n_elem;
    if (!nu || !nv || x.n_elem != nu || y.n_elem != nv ||
            !x.is_finite() || !y.is_finite() || !std::isfinite(cross) ||
            a.has_nan() || b.has_nan() || arma::any(a == R_PosInf) || arma::any(b == R_PosInf))
        Rcpp::stop("Bilinear reduction inputs have invalid dimensions or nonfinite values");
    double high = R_NegInf;
    for (arma::uword j = 0; j < nv; ++j) {
        if (!(j % 64)) Rcpp::checkUserInterrupt();
        const double cy = cross*y[j];
        for (arma::uword i = 0; i < nu; ++i) {
            const double term = (a[i]+b[j])+x[i]*cy;
            if (std::isnan(term) || term == R_PosInf)
                Rcpp::stop("Bilinear reduction has a nonfinite log term");
            if (term > high) high = term;
        }
    }
    if (!std::isfinite(high)) Rcpp::stop("Bilinear reduction has no finite mass");
    arma::vec mu(moments ? nu : 0, arma::fill::zeros), mv(moments ? nv : 0, arma::fill::zeros);
    arma::vec vy(moments ? nu : 0, arma::fill::zeros), ux(moments ? nv : 0, arma::fill::zeros);
    arma::vec cmu = mu, cmv = mv, cvy = vy, cux = ux;
    double total = 0, correction = 0;
    for (arma::uword j = 0; j < nv; ++j) {
        if (!(j % 64)) Rcpp::checkUserInterrupt();
        const double cy = cross*y[j];
        for (arma::uword i = 0; i < nu; ++i) {
            const double weight = std::exp(((a[i]+b[j])+x[i]*cy)-high);
            compensated_add(weight, total, correction);
            if (moments) {
                compensated_add(weight, mu[i], cmu[i]);
                compensated_add(weight, mv[j], cmv[j]);
                compensated_add(weight*y[j], vy[i], cvy[i]);
                compensated_add(weight*x[i], ux[j], cux[j]);
            }
        }
    }
    if (!(total > 0) || !std::isfinite(total)) Rcpp::stop("Bilinear reduction mass is invalid");
    if (moments) {mu /= total; mv /= total; vy /= total; ux /= total;}
    return Rcpp::List::create(Rcpp::Named("value") = high+std::log(total),
        Rcpp::Named("mu") = Rcpp::NumericVector(mu.begin(), mu.end()),
        Rcpp::Named("mv") = Rcpp::NumericVector(mv.begin(), mv.end()),
        Rcpp::Named("weighted_v") = Rcpp::NumericVector(vy.begin(), vy.end()),
        Rcpp::Named("weighted_u") = Rcpp::NumericVector(ux.begin(), ux.end()));
}

// Cell-fast log sums and normalized first scores for the same generic bilinear
// terms. R supplies score_a[i,]+score_b[j,]+cross_score*x[i]*y[j]. One-based
// cell labels need not be sorted; u cells vary fastest in the returned matrix.
// Two passes retain O(ncells_a*ncells_b*(p+1)) memory, never node-pair buffers.
// [[Rcpp::export(rng = false)]]
arma::mat bilinear_cells_cpp(const arma::vec& a, const arma::vec& b,
        const arma::vec& x, const arma::vec& y, double cross,
        const arma::mat& score_a, const arma::mat& score_b,
        const arma::vec& cross_score, const arma::vec& cell_a,
        const arma::vec& cell_b, double ncells_a, double ncells_b) {
    const arma::uword nu = a.n_elem, nv = b.n_elem, p = score_a.n_cols;
    if (!std::isfinite(ncells_a) || !std::isfinite(ncells_b) ||
            ncells_a < 0 || ncells_b < 0 || ncells_a != std::floor(ncells_a) ||
            ncells_b != std::floor(ncells_b) || ncells_a > INT_MAX ||
            ncells_b > INT_MAX || ncells_a*ncells_b > INT_MAX ||
            x.n_elem != nu || y.n_elem != nv || cell_a.n_elem != nu ||
            cell_b.n_elem != nv || score_a.n_rows != nu || score_b.n_rows != nv ||
            score_b.n_cols != p || cross_score.n_elem != p)
        Rcpp::stop("Bilinear cell reduction dimensions do not agree");
    if (!x.is_finite() || !y.is_finite() || !std::isfinite(cross) ||
            a.has_nan() || b.has_nan() || arma::any(a == R_PosInf) ||
            arma::any(b == R_PosInf))
        Rcpp::stop("Bilinear cell reduction inputs have nonfinite values");
    if (!cell_a.is_finite() || !cell_b.is_finite() ||
            arma::any(cell_a != arma::floor(cell_a)) ||
            arma::any(cell_b != arma::floor(cell_b)) ||
            arma::any(cell_a < 1) || arma::any(cell_a > ncells_a) ||
            arma::any(cell_b < 1) || arma::any(cell_b > ncells_b))
        Rcpp::stop("Bilinear cell labels must be one-based integers within their cell counts");
    const arma::uword nca = ncells_a, cells = nca*static_cast<arma::uword>(ncells_b);
    const arma::uvec ca = arma::conv_to<arma::uvec>::from(cell_a-1);
    const arma::uvec cb = arma::conv_to<arma::uvec>::from(cell_b-1);
    arma::vec high(cells); high.fill(R_NegInf);
    arma::mat out(cells, p+1, arma::fill::zeros), correction(cells, p+1, arma::fill::zeros);
    for (arma::uword j = 0; j < nv; ++j) {
        if (!(j % 64)) Rcpp::checkUserInterrupt();
        if (b[j] == R_NegInf) continue;
        const double cy = cross*y[j];
        const arma::uword base = nca*cb[j];
        for (arma::uword i = 0; i < nu; ++i) {
            if (a[i] == R_NegInf) continue;
            const double term = (a[i]+b[j])+x[i]*cy;
            if (std::isnan(term) || term == R_PosInf)
                Rcpp::stop("Bilinear cell reduction has a nonfinite log term");
            const arma::uword cell = base+ca[i];
            if (term > high[cell]) high[cell] = term;
        }
    }
    for (arma::uword j = 0; j < nv; ++j) {
        if (!(j % 64)) Rcpp::checkUserInterrupt();
        if (b[j] == R_NegInf) continue;
        const double cy = cross*y[j];
        const arma::uword base = nca*cb[j];
        for (arma::uword i = 0; i < nu; ++i) {
            if (a[i] == R_NegInf) continue;
            const arma::uword cell = base+ca[i];
            if (high[cell] == R_NegInf) continue;
            const double weight = std::exp(((a[i]+b[j])+x[i]*cy)-high[cell]);
            // Embedded zero weights and numerical underflow mask derivatives.
            if (!weight) continue;
            compensated_add(weight, out(cell, 0), correction(cell, 0));
            for (arma::uword k = 0; k < p; ++k) {
                const double score = (score_a(i,k)+score_b(j,k))+
                    (cross_score[k]*x[i])*y[j];
                if (!std::isfinite(score))
                    Rcpp::stop("Bilinear cell reduction has a nonfinite positive-weight score");
                compensated_add(weight*score, out(cell,k+1), correction(cell,k+1));
            }
        }
    }
    for (arma::uword cell = 0; cell < cells; ++cell) {
        const double mass = out(cell,0);
        if (high[cell] == R_NegInf) {out(cell,0) = R_NegInf; continue;}
        if (!(mass > 0) || !std::isfinite(mass))
            Rcpp::stop("Bilinear cell reduction mass is invalid");
        out(cell,0) = high[cell]+std::log(mass);
        for (arma::uword k = 1; k <= p; ++k) {
            out(cell,k) /= mass;
            if (!std::isfinite(out(cell,k)))
                Rcpp::stop("Bilinear cell reduction normalized score is invalid");
        }
    }
    return out;
}

// Flattened H/T columns follow R's (index, coordinate)/(index, index, coordinate) axes.
// This helper returns centered numerical moments; R assembles log derivatives.
// [[Rcpp::export(rng = false)]]
Rcpp::List mixture_moments_cpp(const arma::vec& terms, const arma::mat& gradient,
        const arma::mat& hessian, const arma::mat& third, double rows, int q_input) {
    const arma::uword d = gradient.n_cols;
    if (q_input < 0) Rcpp::stop("Mixture derivative dimensions do not agree");
    const arma::uword q = q_input;
    if (q > d || (hessian.n_cols && hessian.n_cols != q*d) ||
        (third.n_cols && (third.n_cols != q*q*d || !hessian.n_cols)))
        Rcpp::stop("Mixture derivative dimensions do not agree");
    arma::vec value;
    const arma::mat weights = normalized_weights(terms, rows, value);
    const arma::mat mg = weighted_rows(gradient, weights), mh = weighted_rows(hessian, weights);
    const arma::mat mt = weighted_rows(third, weights);
    const arma::uword n = weights.n_rows;
    arma::mat cg(n, mh.n_cols, arma::fill::zeros), ch(n, mt.n_cols, arma::fill::zeros);
    arma::mat cm = ch, ct = ch;
    for (arma::uword i = 0; i < n && q && mh.n_cols && weights.n_cols; ++i) {
        const arma::vec all = weights.row(i).t();
        const arma::uvec keep = arma::find(all != 0), index = i+n*keep;
        if (keep.is_empty()) continue;
        const arma::vec w = all.elem(keep);
        arma::mat g = gradient.rows(index), h = hessian.rows(index);
        g.each_row() -= mg.row(i);
        h.each_row() -= mh.row(i);
        const arma::mat wg = g.each_col()%w;
        cg.row(i) = arma::vectorise(g.cols(0, q-1).t()*wg).t();
        for (arma::uword j = 0; j < q && mt.n_cols; ++j) for (arma::uword l = 0; l < q; ++l) {
            const arma::uvec columns = arma::regspace<arma::uvec>(j+q*l, q*q, q*q*d-1);
            ch.submat(arma::uvec{i}, columns) = h.col(j+q*l).t()*wg;
            cm.submat(arma::uvec{i}, columns) = (w%g.col(l)).t()*h.cols(arma::regspace<arma::uvec>(j, q, q*d-1));
            ct.submat(arma::uvec{i}, columns) = (g.col(j)%g.col(l)).t()*wg;
        }
    }
    return Rcpp::List::create(Rcpp::Named("value")=Rcpp::NumericVector(value.begin(), value.end()),
        Rcpp::Named("mean_gradient")=mg,
        Rcpp::Named("mean_hessian")=mh, Rcpp::Named("cov_gradient")=cg,
        Rcpp::Named("mean_third")=mt, Rcpp::Named("cov_latent_hessian_score")=ch,
        Rcpp::Named("cov_mixed_hessian_score")=cm, Rcpp::Named("central_third")=ct);
}

// Row contraction: score is n x q and derivative n x (q*p), q varying fastest
// within each of the p blocks; column k sums score[, j]*derivative[, block k, j].
// [[Rcpp::export(rng = false)]]
arma::mat row_contract_cpp(const arma::mat& score, const arma::mat& derivative, int p_input) {
    if (p_input < 0) Rcpp::stop("Contraction dimensions do not agree");
    const arma::uword q = score.n_cols, p = p_input;
    if (derivative.n_rows != score.n_rows || derivative.n_cols != q*p)
        Rcpp::stop("Contraction dimensions do not agree");
    arma::mat out(score.n_rows, p, arma::fill::zeros);
    for (arma::uword k = 0; k < p && q; ++k)
        out.col(k) = arma::sum(score%derivative.cols(k*q, (k+1)*q-1), 1);
    return out;
}
