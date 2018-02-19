// [[Rcpp::depends(RcppArmadillo)]]

# include "ncpen.h"

//' @title
//' Native object function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param name .
//' @param y_vec .
//' @param x_mat .
//' @param b_vec .
//' @param r_eff .
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
double native_cpp_obj_fun_(std::string name, arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
     obj_fun_ptr obj_fun = get_obj_fun_ptr(name);

     return obj_fun(y_vec, x_mat, b_vec, r_eff);
}

//' @title
//' Native object gradient function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param name .
//' @param y_vec .
//' @param x_mat .
//' @param b_vec .
//' @param r_eff .
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
arma::vec native_cpp_obj_grad_fun_(std::string name, arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
     obj_grad_fun_ptr obj_grad_fun = get_obj_grad_fun_ptr(name);

     return obj_grad_fun(y_vec, x_mat, b_vec, r_eff);
}


//' @title
//' Native QLASSO function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param q_mat .
//' @param l_vec .
//' @param b_vec0 .
//' @param w_vec .
//' @param lam .
//' @param iter_max .
//' @param b_eps .
//' @param k_eps .
//' @param p_eff .
//' @param q_rank .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
Rcpp::List native_cpp_qlasso_fun_(arma::mat& q_mat, arma::vec& l_vec, arma::vec& b_vec0, arma::vec& w_vec,
                                   double lam, double iter_max, double b_eps, double k_eps, arma::uword p_eff, arma::uword q_rank) {
     p_ncpen_ret ret_buff;
     qlasso_fun(q_mat, l_vec, b_vec0, w_vec,lam, iter_max, b_eps, k_eps, p_eff, q_rank, ret_buff);

     return Rcpp::List::create(Rcpp::Named("g.vec") = ret_buff.g_vec,
                               Rcpp::Named("b.vec") = ret_buff.b_vec,
                               Rcpp::Named("f.vec") = ret_buff.f_vec,
                               Rcpp::Named("con") = ret_buff.con);


}

//' @title
//' Native point ncpen function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param y_vec .
//' @param x_mat .
//' @param b_vec .
//' @param w_vec .
//' @param lam .
//' @param gam .
//' @param tau .
//' @param iter_max .
//' @param b_eps .
//' @param k_eps .
//' @param p_eff .
//' @param r_eff .
//' @param family .
//' @param penalty .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
Rcpp::List native_cpp_p_ncpen_fun_(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, arma::vec& w_vec,
                                   double lam, double gam, double tau, double iter_max, double b_eps, double k_eps, arma::uword p_eff, double r_eff,
                                   SEXP family, SEXP penalty) {
     // Rcpp::List p_ncpen_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec0, arma::vec& w_vec,
     //                        double lam, double gam, double tau, double iter_max, double b_eps, double k_eps, arma::uword p_eff, double r_eff,
     //                        pen_fun_ptr pen_fun, pen_grad_fun_ptr pen_grad_fun,
     //                        obj_fun_ptr obj_fun, obj_grad_fun_ptr obj_grad_fun, obj_hess_fun_ptr obj_hess_fun) {

     std::string fam = Rcpp::as<std::string>(family);
     std::string pen = Rcpp::as<std::string>(penalty);

     // return p_ncpen_fun(y_vec, x_mat, b_vec, w_vec, lam, gam, tau, iter_max, b_eps, k_eps, p_eff, r_eff,
     //                    get_obj_fun_ptr(fam), get_obj_grad_fun_ptr(fam), get_obj_hess_fun_ptr(fam),
     //                    get_pen_fun_ptr(pen), get_pen_grad_fun_ptr(pen));

     p_ncpen_ret ret_buff;
     p_ncpen_fun(y_vec, x_mat, b_vec, w_vec, lam, gam, tau, iter_max, b_eps, k_eps, p_eff, r_eff,
                 get_obj_fun_ptr(fam), get_obj_grad_fun_ptr(fam), get_obj_hess_fun_ptr(fam),
                 get_pen_fun_ptr(pen), get_pen_grad_fun_ptr(pen), ret_buff);

     return Rcpp::List::create(Rcpp::Named("g.vec") = ret_buff.g_vec,
                               Rcpp::Named("b.vec") = ret_buff.b_vec,
                               Rcpp::Named("f.vec") = ret_buff.f_vec,
                               Rcpp::Named("con") = ret_buff.con);

}

//' @title
//' Native point ncpen function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param y_vec .
//' @param x_mat0 .
//' @param x_std .
//' @param intc .
//' @param w_vec0 .
//' @param lam_vec0 .
//' @param r_lam .
//' @param gam .
//' @param tau .
//' @param p_max .
//' @param iter_max .
//' @param b_eps .
//' @param k_eps .
//' @param p_eff .
//' @param r_eff .
//' @param family .
//' @param penalty .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
Rcpp::List native_cpp_ncpen_fun_(arma::vec& y_vec, arma::mat& x_mat0, bool x_std, bool intc,
                                 arma::vec& w_vec0, arma::vec& lam_vec0, double r_lam,
                                 double gam, double tau, arma::uword p_max, double iter_max, double b_eps, double k_eps,
                                 arma::uword p_eff, double r_eff,
                                 SEXP family, SEXP penalty) {

     std::string fam = Rcpp::as<std::string>(family);
     std::string pen = Rcpp::as<std::string>(penalty);

     ncpen_ret ret_buff;
     ncpen_fun(y_vec, x_mat0, x_std, intc,
               w_vec0, lam_vec0, r_lam,
               gam, tau, p_max, iter_max, b_eps, k_eps,
               p_eff, r_eff,
               fam, pen, ret_buff);

     return Rcpp::List::create(Rcpp::Named("b.mat") = ret_buff.b_mat,
                               Rcpp::Named("g.mat") = ret_buff.g_mat,
                               Rcpp::Named("f.vec") = ret_buff.f_vec,
                               Rcpp::Named("c.mat") = ret_buff.c_mat,
                               Rcpp::Named("lam.vec") = ret_buff.lam_vec,
                               Rcpp::Named("d.vec") = ret_buff.d_vec,
                               Rcpp::Named("w.vec") = ret_buff.w_vec,
                               Rcpp::Named("warnings") = ret_buff.warnings);

}

