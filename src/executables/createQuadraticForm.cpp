//
// Created by Lars Gebraad on 11/06/18.
//

#include <cstdlib>
#include <armadillo>
#include "../hmc/linearSampler.hpp"

double X(arma::mat A, arma::colvec B, arma::mat C, arma::vec m) {
    return arma::as_scalar(m.t() * A * m + B.t() * m + C);
}

int main(int argc, char *argv[]) {

    arma::mat G;
    arma::mat d0;

    G.load(argv[1]);
    arma::vec m0 = 6.67e-4 * ones(G.n_cols, 1);
    d0.load(argv[2]);
    m0.save(argv[3], raw_ascii);

    std::stringstream m_var_s(argv[4]);
    double m_var;
    m_var_s >> m_var;
    std::stringstream d_var_s(argv[5]);
    double d_var;
    d_var_s >> d_var;

    arma::sp_mat G_sp = sp_mat(G);
    arma::sp_mat cm = m_var * speye<sp_mat>(m0.n_rows, m0.n_rows);
    arma::sp_mat cd = d_var * speye<sp_mat>(d0.n_rows, d0.n_rows);

    arma::mat A;
    arma::colvec B;
    arma::mat C;

    arma::sp_mat invcm = speye<sp_mat>(m0.n_rows, m0.n_rows) / (m_var);
    arma::sp_mat invcd = speye<sp_mat>(d0.n_rows, d0.n_rows) / (d_var);

    arma::mat invcdG = invcd * G;

    A = 0.5 * (invcm + G.t() * invcdG);
    B = -(m0.t() * invcm + d0.t() * invcdG).t();
    C = 0.5 * (m0.t() * invcm * m0 + d0.t() * invcd * d0);

    A.save("A.txt", raw_ascii);
    B.save("B.txt", raw_ascii);
    C.save("C.txt", raw_ascii);

    arma::vec m_post = m0 + (cm * G.t()) * inv(G * cm * G.t() + cd) * (d0 - G * m0);

    arma::vec m_post_alt = -arma::inv(A) * B / 2;

//    cout << m_post - m_post_alt;

    arma::mat cm_post = inv(G.t() * invcd * G + invcm);

//    cout << cm_post - inv(2 * A);

    cout << X(A, B, C, m_post);

    m_post.save("m_post.txt", raw_ascii);
    cm_post.save("cm_post.txt", raw_ascii);

    return EXIT_SUCCESS;
}
