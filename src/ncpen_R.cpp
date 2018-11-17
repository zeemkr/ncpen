// [[Rcpp::depends(RcppArmadillo)]]

# include "ncpen.h"

//' @title
//' N/A.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param dev_mode .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
int native_cpp_set_dev_mode_(bool dev_mode) {
     return set_dev_mode(dev_mode);
}

//' @title
//' N/A.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param fam .
//' @param y_vec .
//' @param x_mat .
//' @param iter_max .
//' @param b_eps .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
arma::vec native_cpp_nr_fun_(std::string fam, arma::vec& y_vec, arma::mat& x_mat, double iter_max, double b_eps) {
     return nr_fun(fam, y_vec, x_mat, iter_max, b_eps);
}

//' @title
//' Native Penalty function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param name .
//' @param b_vec .
//' @param lam .
//' @param gam .
//' @param tau .
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
arma::vec native_cpp_pen_fun_(std::string name, arma::vec& b_vec, double lam, double gam, double tau) {
     pen_fun_ptr pen_fun = get_pen_fun_ptr(name);

     return pen_fun(b_vec, lam, gam, tau);
}

//' @title
//' Native Penalty Gradient function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param name .
//' @param b_vec .
//' @param lam .
//' @param gam .
//' @param tau .
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
arma::vec native_cpp_pen_grad_fun_(std::string name, arma::vec& b_vec, double lam, double gam, double tau) {
     pen_grad_fun_ptr pen_grad_fun = get_pen_grad_fun_ptr(name);

     return pen_grad_fun(b_vec, lam, gam, tau);
}


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
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
double native_cpp_obj_fun_(std::string name, arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     obj_fun_ptr obj_fun = get_obj_fun_ptr(name);

     return obj_fun(y_vec, x_mat, b_vec);
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
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
arma::vec native_cpp_obj_grad_fun_(std::string name, arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     obj_grad_fun_ptr obj_grad_fun = get_obj_grad_fun_ptr(name);

     return obj_grad_fun(y_vec, x_mat, b_vec);
}

//' @title
//' Native object Hessian function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param name .
//' @param y_vec .
//' @param x_mat .
//' @param b_vec .
//'
//'
//' @return
//' .
//'
// [[Rcpp::export]]
arma::mat native_cpp_obj_hess_fun_(std::string name, arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     obj_hess_fun_ptr obj_hess_fun = get_obj_hess_fun_ptr(name);

     return obj_hess_fun(y_vec, x_mat, b_vec);
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
//' @param iiter_max .
//' @param b_eps .
//' @param k_eps .
//' @param p_eff .
//' @param q_rank .
//' @param cut .
//' @param c_eps .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
Rcpp::List native_cpp_qlasso_fun_(arma::mat& q_mat, arma::vec& l_vec, arma::vec& b_vec0, arma::vec& w_vec,
                                  double lam, double iter_max, double iiter_max, double b_eps, double k_eps,
                                  arma::uword p_eff, arma::uword q_rank, bool cut, double c_eps) {
     p_ncpen_ret ret_buff;
     qlasso_fun(q_mat, l_vec, b_vec0, w_vec,lam, iter_max, iiter_max, b_eps, k_eps, p_eff, q_rank, cut, c_eps, ret_buff);

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
//' @param alp .
//' @param iter_max .
//' @param qiter_max .
//' @param qiiter_max .
//' @param b_eps .
//' @param k_eps .
//' @param p_eff .
//' @param cut .
//' @param c_eps .
//' @param family .
//' @param penalty .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
Rcpp::List native_cpp_p_ncpen_fun_(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, arma::vec& w_vec,
                                   double lam, double gam, double tau, double alp, double iter_max, double qiter_max, double qiiter_max,
                                   double b_eps, double k_eps, arma::uword p_eff, bool cut, double c_eps,
                                   SEXP family, SEXP penalty) {
     std::string fam = Rcpp::as<std::string>(family);
     std::string pen = Rcpp::as<std::string>(penalty);

     p_ncpen_ret ret_buff;
     p_ncpen_fun(y_vec, x_mat, b_vec, w_vec, lam, gam, tau, alp, iter_max, qiter_max,
                 qiiter_max, b_eps, k_eps, p_eff, cut, c_eps,
                 get_obj_fun_ptr(fam), get_obj_grad_fun_ptr(fam), get_obj_hess_fun_ptr(fam),
                 get_pen_fun_ptr(pen), get_pen_grad_fun_ptr(pen), ret_buff);

     return Rcpp::List::create(Rcpp::Named("g.vec") = ret_buff.g_vec,
                               Rcpp::Named("b.vec") = ret_buff.b_vec,
                               Rcpp::Named("f.vec") = ret_buff.f_vec,
                               Rcpp::Named("con") = ret_buff.con);

}

//' @title
//' Native ncpen function.
//'
//' @description
//' This is internal use only function. Manual left blank on purpose.
//'
//' @param y_vec .
//' @param x_mat0 .
//' @param w_vec0 .
//' @param lam_vec0 .
//' @param gam .
//' @param tau .
//' @param alp .
//' @param d_max .
//' @param iter_max .
//' @param qiter_max .
//' @param qiiter_max .
//' @param b_eps .
//' @param k_eps .
//' @param p_eff .
//' @param cut .
//' @param c_eps .
//' @param add .
//' @param family .
//' @param penalty .
//' @param loc .
//' @param ob_vec .
//' @param div .
//'
//' @return
//' .
//'
// [[Rcpp::export]]
Rcpp::List native_cpp_ncpen_fun_(arma::vec& y_vec, arma::mat& x_mat0,
                                 arma::vec& w_vec0, arma::vec& lam_vec0,
                                 double gam, double tau, double alp, arma::uword d_max, double iter_max, double qiter_max, double qiiter_max,
                                 double b_eps, double k_eps,
                                 arma::uword p_eff, bool cut, double c_eps, arma::uword add,
                                 SEXP family, SEXP penalty,
                                 bool loc, arma::vec& ob_vec, int div) {

     std::string fam = Rcpp::as<std::string>(family);
     std::string pen = Rcpp::as<std::string>(penalty);

     ncpen_ret ret_buff;
     ncpen_fun(y_vec, x_mat0,w_vec0, lam_vec0,
               gam, tau, alp, d_max, iter_max, qiter_max, qiiter_max, b_eps, k_eps,
               p_eff, cut, c_eps, add,
               fam, pen, loc, ob_vec, div, ret_buff);

     return Rcpp::List::create(Rcpp::Named("beta") = ret_buff.b_mat,
                               Rcpp::Named("grad") = ret_buff.g_mat,
                               Rcpp::Named("f.vec") = ret_buff.f_vec,
                               Rcpp::Named("conv") = ret_buff.c_mat,
                               Rcpp::Named("lambda") = ret_buff.lam_vec,
                               Rcpp::Named("df") = ret_buff.d_vec,
                               Rcpp::Named("w.lambda") = ret_buff.w_vec,
                               Rcpp::Named("warnings") = ret_buff.warnings);

}

