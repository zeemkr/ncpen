#include <RcppArmadillo.h>

/*--------------------------------------------------------
 * Type definitions
 --------------------------------------------------------*/
const int CONST_N_BIG_EXP_LOADER = 700;
const double POI_BIG_XB = 700;
const double LOGI_BIG_XB = 700;
const double LOGI_SMALL_NUMBER = 1e-7;

/* End of type definitions --------------------------------------------------------*/

// qlasso and p_ncpen return buffer
class p_ncpen_ret {
public:
     arma::vec g_vec;
     arma::vec b_vec;
     arma::vec f_vec;
     bool con;
};

// ncpen return buffer
class ncpen_ret {
public:
     arma::mat b_mat;
     arma::mat g_mat;
     arma::vec f_vec;
     arma::mat c_mat;
     arma::vec lam_vec;
     arma::vec d_vec;
     arma::vec w_vec;
     arma::uvec warnings;
};

// Loss function proto types
typedef double (*obj_fun_ptr)(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff);
typedef arma::vec (*obj_grad_fun_ptr)(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff);
typedef arma::mat (*obj_hess_fun_ptr)(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff);

obj_fun_ptr get_obj_fun_ptr(std::string name);
obj_grad_fun_ptr get_obj_grad_fun_ptr(std::string name);
obj_hess_fun_ptr get_obj_hess_fun_ptr(std::string name);

// Penalty function proto type
typedef arma::vec (*pen_fun_ptr)(arma::vec& b_vec, double lam, double gam, double tau);
typedef arma::vec (*pen_grad_fun_ptr)(arma::vec& b_vec, double lam, double gam, double tau);

pen_fun_ptr get_pen_fun_ptr(std::string name);
pen_grad_fun_ptr get_pen_grad_fun_ptr(std::string name);

int qlasso_fun(arma::mat& q_mat, arma::vec& l_vec, arma::vec& b_vec0, arma::vec& w_vec,
               double lam, double iter_max, double b_eps, double k_eps, arma::uword p_eff, arma::uword q_rank, p_ncpen_ret& ret_buff);

int p_ncpen_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec0, arma::vec& w_vec,
                double lam, double gam, double tau, double iter_max, double b_eps, double k_eps, arma::uword p_eff, double r_eff,
                obj_fun_ptr obj_fun, obj_grad_fun_ptr obj_grad_fun, obj_hess_fun_ptr obj_hess_fun,
                pen_fun_ptr pen_fun, pen_grad_fun_ptr pen_grad_fun, p_ncpen_ret& ret_buff);

int ncpen_fun(arma::vec& y_vec, arma::mat& x_mat0, bool x_std, bool intc,
              arma::vec& w_vec0, arma::vec& lam_vec0, double r_lam,
              double gam, double tau, arma::uword p_max, double iter_max, double b_eps, double k_eps,
              arma::uword p_eff, double r_eff,
              std::string fam, std::string pen, ncpen_ret& ret_buff);


