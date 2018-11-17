/*--------------------------------------------------------
 * Includes
 --------------------------------------------------------*/
#include "ncpen.h"
#include <RcppArmadillo.h>

using namespace arma;

/* End of includes --------------------------------------------------------*/

// Global functions ----------------------------
bool NCPEN_DEVELOP_MODE = false;

int set_dev_mode(bool dev_mode) {
     NCPEN_DEVELOP_MODE = dev_mode;

     return 0;
}
//----------------------------------------------------

/*--------------------------------------------------------
* Utility functions
--------------------------------------------------------*/
// This is for a scalar number.
template <typename T> int sign_scalar_(T val) {
     return (T(0) < val) - (val < T(0));
}

// arma::vec row_pmin(const arma::vec& v1, const double d1) {
// //     arma::mat pmin_mat(v1.n_rows, 2);
// //     pmin_mat.col(0) = v1;
// //     pmin_mat.col(1).fill(d1);
// //     return arma::min(pmin_mat, 1); // for each row
//
//      arma::vec v = v1;
//      v(arma::find(v>d1)).fill(d1);
//      return v;
// }

// arma::vec rm_row(arma::uword rm_idx, arma::vec& v) {
//      //     Rcout << "rm_idx: " << rm_idx << std::endl;
//      arma::uword len = v.n_rows;
//      arma::vec ret(len-1);
//      ret.head(rm_idx) = v.head(rm_idx);
//      ret.tail(len - rm_idx - 1) = v.tail(len - rm_idx - 1);
//
//      return ret;
// }


// arma::vec rm_row(arma::vec rm_vec, arma::uword start, arma::uword end) {
//      rm_vec = sort(rm_vec);
//      arma::uword len = (end-start + 1) - rm_vec.n_rows;
//      arma::vec ret(len);
//
//      arma::uword ii = 0;
//      arma::uword rm_i = 0;
//      for(arma::uword i = 0; i<len; i++) {
//           if(rm_vec(rm_i) != i) {
//                ret(ii) = i;
//                ii++;
//           } else {
//                rm_i++;
//           }
//      }
//
//      return ret;
// }

// int which_max(arma::vec& v) { //do not use uword to return -1
//      if(v.n_rows <=0) return -1;
//
//      arma::uword max_idx = 0;
//      double max_val = v[0];
//
//      for(arma::uword i=1; i<v.n_rows; i++) {
//           if(v(i) > max_val) {
//                max_val = v(i);
//                max_idx = i;
//           }
//      }
//
//      return max_idx;
// }

// arma::uvec rm_row(arma::uvec rm_vec, arma::uword start, arma::uword end) {
//      rm_vec = sort(rm_vec);
//      arma::uword len = (end-start + 1) - rm_vec.n_rows;
//      arma::uvec ret(len);
//
//      arma::uword ii = 0;
//      arma::uword rm_i = 0;
//      for(arma::uword i = 0; ; i++) {
//           //Rcout << "i: " <<i<<"  ii: " <<ii<<"  rm_i: " <<rm_i <<endl;
//           if(rm_i == rm_vec.n_rows) {
//                for(;ii<len;) {
//                     ret(ii) = i;
//                     ii++;
//                     i++;
//                }
//                break;
//           }
//
//           if(rm_vec(rm_i) != i) {
//                ret(ii) = i;
//                ii++;
//           } else {
//                rm_i++;
//           }
//      }
//
//      return ret;
// }

// call by value
// arma::vec arma_u2d(arma::uvec& v) {
//      return arma::ones<vec>(v.n_rows)%v;
// }

// arma::uvec arma_or(arma::uvec& a, arma::uvec& b) {
//      return (a + b) > 0;
// }
//
// arma::uvec arma_and(arma::uvec& a, arma::uvec& b) {
//      return (a % b) > 0;
// }


/* End of utility functions --------------------------------------------------------*/



/*--------------------------------------------------------
* Main functions
--------------------------------------------------------*/

// 1. Penalty functions ------------------------------------------

// if(penalty=="lasso"){
//      pen.fun = function(b.vec,lam,gam,tau){
arma::vec lasso_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     //           return(lam*abs(b.vec))
     return lam*arma::abs(b_vec);
}
//      pen.grad.fun = function(b.vec,lam,gam,tau){
arma::vec lasso_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     //           return(lam*sign(b.vec))
     return lam*arma::sign(b_vec);
}

// if(pen=="scad"){
//      pen.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec)
//           tem0.vec = (lam*ab.vec) * (ab.vec<lam)
//           tem1.vec = ((tau*lam*(ab.vec-lam)-(ab.vec^2-lam^2)/2)/(tau-1)+lam^2) * (ab.vec>=lam)*(ab.vec<tau*lam)
//           tem2.vec = ((tau+1)*lam^2/2) * (ab.vec>=tau*lam)
//           return(tem0.vec+tem1.vec+tem2.vec)
//      }
//      pen.grad.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec); sb.vec = sign(b.vec)
//                tem0.vec = (lam) * (ab.vec<lam)
//                tem1.vec = ((tau*lam-ab.vec)/(tau-1)) * (ab.vec>=lam)*(ab.vec<tau*lam)
//                return((tem0.vec+tem1.vec)*sb.vec)
//      }
// }
arma::vec scad_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec tem0_vec = (lam*ab_vec) % (ab_vec < lam);
     arma::vec tem1_vec = ((tau*lam*(ab_vec-lam)-(arma::pow(ab_vec,2)-std::pow(lam,2))/2)/(tau-1)+std::pow(lam,2)) % (ab_vec>=lam)%(ab_vec<tau*lam);
     //     arma::vec tem2_vec = ((tau+1)*std::pow(lam, 2)/2) * arma_u2d(ab_vec>=tau*lam);
     arma::vec tem2_vec = ((tau+1)*std::pow(lam, 2)/2) * conv_to<vec>::from(ab_vec>=tau*lam);

     return tem0_vec + tem1_vec + tem2_vec;
}

arma::vec scad_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     //     arma::vec sb_vec = arma::sign(b_vec);
     arma::vec tem0_vec = lam * conv_to<vec>::from(ab_vec<lam);
     arma::vec tem1_vec = ((tau*lam-ab_vec)/(tau-1)) % (ab_vec>=lam)%(ab_vec<tau*lam);

     //     return (tem0_vec+tem1_vec)%sb.vec;
     return (tem0_vec+tem1_vec)%arma::sign(b_vec);
}

// if(pen=="mcp"){
//      pen.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec)
//           tem0.vec = (lam*ab.vec-ab.vec^2/2/tau) * (ab.vec<tau*lam)
//           tem1.vec = (tau*lam^2/2) * (ab.vec>=tau*lam)
//           return(tem0.vec+tem1.vec)
//      }
//      pen.grad.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec); sb.vec = sign(b.vec)
//                tem0.vec = (lam-ab.vec/tau) * (ab.vec<tau*lam)
//                return(tem0.vec*sb.vec)
//      }
// }
arma::vec mcp_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec tem0_vec = (lam*ab_vec-arma::pow(ab_vec,2)/2/tau) % (ab_vec<tau*lam);
     arma::vec tem1_vec = (tau*std::pow(lam, 2)/2) * conv_to<vec>::from(ab_vec>=tau*lam);

     return tem0_vec + tem1_vec;
}

arma::vec mcp_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     //     arma::vec sb_vec = arma::sign(b_vec);
     arma::vec tem0_vec = (lam-ab_vec/tau) % (ab_vec<tau*lam);

     //     return tem0_vec%sb_vec;
     return tem0_vec%arma::sign(b_vec);
}

// if(pen=="tlp"){
//      pen.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec)
//           tem0.vec = (lam*ab.vec) * (ab.vec<tau)
//           tem1.vec = (lam*tau) * (ab.vec>=tau)
//           return(tem0.vec+tem1.vec)
//      }
//      pen.grad.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec); sb.vec = sign(b.vec)
//                tem0.vec = (lam) * (ab.vec<tau)
//                return(tem0.vec*sb.vec)
//      }
// }

arma::vec tlp_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec tem0_vec = (lam*ab_vec) % (ab_vec<tau);
     arma::vec tem1_vec = (lam*tau) * conv_to<vec>::from(ab_vec>=tau);

     return tem0_vec + tem1_vec;
}

arma::vec tlp_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     //     arma::vec sb_vec = arma::sign(b_vec);
     arma::vec tem0_vec = (lam) * conv_to<vec>::from(ab_vec<tau);

     //     return tem0_vec%sb.vec;
     return tem0_vec%arma::sign(b_vec);
}
//
// if(pen=="classo"){
//      pen.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec)
//           tem0.vec = (-ab.vec^2/tau/2+lam*ab.vec) * (ab.vec<tau*(lam-gam))
//           tem1.vec = (gam*ab.vec-tau^2*(lam-gam)^2/tau/2+lam*tau*(lam-gam)-tau*gam*(lam-gam)) * (ab.vec>=tau*(lam-gam))
//           return(tem0.vec+tem1.vec)
//      }
//      pen.grad.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec); sb.vec = sign(b.vec)
//                tem0.vec = (lam-ab.vec/tau) * (ab.vec<tau*(lam-gam))
//                tem1.vec = (gam) * (ab.vec>=tau*(lam-gam))
//                return((tem0.vec+tem1.vec)*sb.vec)
//      }
// }
arma::vec classo_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec tem0_vec = (-arma::pow(ab_vec, 2)/tau/2+lam*ab_vec) % (ab_vec<tau*(lam-gam));
     arma::vec tem1_vec = (gam*ab_vec-std::pow(tau,2)*std::pow(lam-gam,2)/tau/2+lam*tau*(lam-gam)-tau*gam*(lam-gam)) % (ab_vec>=tau*(lam-gam));

     return tem0_vec + tem1_vec;
}

arma::vec classo_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     //     arma::vec sb_vec = arma::sign(b_vec);
     arma::vec tem0_vec = (lam-ab_vec/tau) % (ab_vec<tau*(lam-gam));
     arma::vec tem1_vec = (gam) * conv_to<vec>::from(ab_vec>=tau*(lam-gam));

     //     return (tem0_vec+tem1_vec)%sb_vec;
     return (tem0_vec+tem1_vec)%arma::sign(b_vec);
}

// if(pen=="sridge"){ ### modified on 2018 0707
//      pen.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec)
//           tem0.vec = (-ab.vec^2/2/tau+lam*ab.vec)*(ab.vec<tau*lam/(1+tau*gam))
//           tem1.vec = (gam*ab.vec^2/2+tau*lam^2/(1+tau*gam)/2)*(ab.vec>=tau*lam/(1+tau*gam))
//           return(tem0.vec+tem1.vec)
//      }
//      pen.grad.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec); sb.vec = sign(b.vec)
//                tem0.vec = (lam-ab.vec/tau)*(ab.vec<tau*lam/(1+tau*gam))
//                tem1.vec = gam*ab.vec*(ab.vec>=tau*lam/(1+tau*gam))
//                return((tem0.vec+tem1.vec)*sb.vec)
//      }
// }

// pen.fun = function(b.vec,lam,gam,tau){
arma::vec sridge_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     // ab.vec = abs(b.vec)
     arma::vec ab_vec = arma::abs(b_vec);
     // tem0.vec = (-ab.vec^2/2/tau+lam*ab.vec)*(ab.vec<tau*lam/(1+tau*gam))
     arma::vec tem0_vec = (-arma::pow(ab_vec,2)/2/tau+lam*ab_vec)%(ab_vec<tau*lam/(1+tau*gam));
     // tem1.vec = (gam*ab.vec^2/2+tau*lam^2/(1+tau*gam)/2)*(ab.vec>=tau*lam/(1+tau*gam))
     arma::vec tem1_vec = (gam*arma::pow(ab_vec,2)/2+tau*std::pow(lam,2)/(1+tau*gam)/2)%(ab_vec>=tau*lam/(1+tau*gam));
     // return(tem0.vec+tem1.vec)
     return tem0_vec+tem1_vec;
}

// pen.grad.fun = function(b.vec,lam,gam,tau){
arma::vec sridge_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     // ab.vec = abs(b.vec); sb.vec = sign(b.vec)
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec sb_vec = arma::sign(b_vec);
     // tem0.vec = (lam-ab.vec/tau)*(ab.vec<tau*lam/(1+tau*gam))
     arma::vec tem0_vec = (lam-ab_vec/tau)%(ab_vec<tau*lam/(1+tau*gam));
     // tem1.vec = gam*ab.vec*(ab.vec>=tau*lam/(1+tau*gam))
     arma::vec tem1_vec = gam*ab_vec%(ab_vec>=tau*lam/(1+tau*gam));
     // return((tem0.vec+tem1.vec)*sb.vec)
     return (tem0_vec+tem1_vec)%sb_vec;
}


// if(pen=="mbridge"){ Add: 8/31/2016
arma::vec mbridge_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     // ab.vec = abs(b.vec)
     arma::vec ab_vec = arma::abs(b_vec);
     // pen.vec = rep(0,length(b.vec))
     arma::vec pen_vec = arma::zeros<arma::vec>(b_vec.n_rows);
     // pen.vec[ab.vec< tau] = lam*ab.vec[ab.vec<tau]
     arma::uvec idx1 = arma::find(ab_vec < tau);
     pen_vec(idx1) = lam*ab_vec(idx1);
     // pen.vec[ab.vec>=tau] = lam*(2*sqrt(tau*ab.vec[ab.vec>=tau])-tau)
     arma::uvec idx2 = arma::find(ab_vec >= tau);
     pen_vec(idx2) = lam*(2*arma::sqrt(tau*ab_vec(idx2)-tau));
     // return(pen.vec)
     return pen_vec;
}

arma::vec mbridge_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     // ab.vec = abs(b.vec); sb.vec = sign(b.vec)
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec sb_vec = arma::sign(b_vec);
     // grad.vec = rep(0,length(b.vec))
     arma::vec grad_vec = arma::zeros<arma::vec>(b_vec.n_rows);
     // grad.vec[ab.vec< tau] = lam
     grad_vec(arma::find(ab_vec < tau)).fill(lam);
     // grad.vec[ab.vec>=tau] = lam*sqrt(tau/ab.vec[ab.vec>=tau])
     arma::uvec idx2 = arma::find(ab_vec >= tau);
     grad_vec(idx2) = lam*arma::sqrt(tau/ab_vec(idx2));
     // return(grad.vec*sb.vec)
     return(grad_vec%sb_vec);

}

// if(pen=="mlog"){ Add: 8/31/2016
arma::vec mlog_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     // ab.vec = abs(b.vec)
     arma::vec ab_vec = arma::abs(b_vec);
     // pen.vec = rep(0,length(b.vec))
     arma::vec pen_vec = arma::zeros<arma::vec>(b_vec.n_rows);
     // pen.vec[ab.vec< tau] = lam*ab.vec[ab.vec<tau]
     arma::uvec idx1 = arma::find(ab_vec < tau);
     pen_vec(idx1) = lam*ab_vec(idx1);
     // pen.vec[ab.vec>=tau] = lam*tau*(1+log(ab.vec[ab.vec>=tau]/tau))
     arma::uvec idx2 = arma::find(ab_vec >= tau);
     pen_vec(idx2) = lam*tau*(1+arma::log(ab_vec(idx2)/tau));
     // return(pen.vec)
     return pen_vec;
}

arma::vec mlog_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     // ab.vec = abs(b.vec); sb.vec = sign(b.vec)
     arma::vec ab_vec = arma::abs(b_vec);
     arma::vec sb_vec = arma::sign(b_vec);
     // grad.vec = rep(0,length(b.vec))
     arma::vec grad_vec = arma::zeros<arma::vec>(b_vec.n_rows);
     // grad.vec[ab.vec< tau] = lam
     grad_vec(arma::find(ab_vec < tau)).fill(lam);
     // grad.vec[ab.vec>=tau] = lam*tau/ab.vec[ab.vec>=tau]
     arma::uvec idx2 = arma::find(ab_vec >= tau);
     grad_vec(idx2) = lam*tau/ab_vec(idx2);
     // return(grad.vec*sb.vec)
     return(grad_vec%sb_vec);
}

pen_fun_ptr get_pen_fun_ptr(std::string name) {
     if(name.compare("scad") == 0) {
          return scad_pen_fun;
     } else if(name.compare("mcp") == 0) {
          return mcp_pen_fun;
     } else if(name.compare("tlp") == 0) {
          return tlp_pen_fun;
     } else if(name.compare("classo") == 0) {
          return classo_pen_fun;
     } else if(name.compare("sridge") == 0) {
          return sridge_pen_fun;
     } else if(name.compare("mbridge") == 0) {
          return mbridge_pen_fun;
     } else if(name.compare("mlog") == 0) {
          return mlog_pen_fun;
     } else if(name.compare("lasso") == 0) {
          return lasso_pen_fun;
     } else if(name.compare("ridge") == 0) {
          return scad_pen_fun;
     } else {
          throw std::invalid_argument("Invalid penalty funtion option. Only available \"scad\", \"mcp\", \"tlp\", \"classo\", \"sridge\", \"mbridge\", \"mlog\" or \"lasso\".");
          return NULL;
     }
}

pen_grad_fun_ptr get_pen_grad_fun_ptr(std::string name) {
     if(name.compare("scad") == 0) {
          return scad_pen_grad_fun;
     } else if(name.compare("mcp") == 0) {
          return mcp_pen_grad_fun;
     } else if(name.compare("tlp") == 0) {
          return tlp_pen_grad_fun;
     } else if(name.compare("classo") == 0) {
          return classo_pen_grad_fun;
     } else if(name.compare("sridge") == 0) {
          return sridge_pen_grad_fun;
     } else if(name.compare("mbridge") == 0) {
          return mbridge_pen_grad_fun;
     } else if(name.compare("mlog") == 0) {
          return mlog_pen_grad_fun;
     } else if(name.compare("lasso") == 0) {
          return lasso_pen_grad_fun;
     } else if(name.compare("ridge") == 0) {
          return scad_pen_grad_fun;
     } else {
          throw std::invalid_argument("Invalid penalty gradient funtion option. Only available \"scad\", \"mcp\", \"tlp\", \"classo\", \"sridge\", \"mbridge\", \"mlog\" or \"lasso\".");
          return NULL;
     }
}
// End of penalty functions ------------------------------------------


// 2. Loss functions ------------------------------------------

// 1. if(fam=="lin"){
double lin_obj_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = drop(x.mat%*%b.vec)
     arma::vec xb_vec = x_mat*b_vec;
     //   return(  sum((xb.vec-y.vec)^2)/length(y.vec)/2)
     return arma::sum(arma::pow(xb_vec-y_vec, 2))/y_vec.n_rows/2; //(Modified) 20180707
}

arma::vec lin_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = drop(x.mat%*%b.vec)
     arma::vec xb_vec = x_mat*b_vec;
     //   return( drop(t(x.mat)%*%(xb.vec-y.vec)/length(y.vec)))
     return x_mat.t()*(xb_vec-y_vec)/y_vec.n_rows;  //(Modified) 20180707
}

//
arma::mat lin_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   return( t(x.mat)%*%x.mat/length(y.vec))
     return x_mat.t()*x_mat/y_vec.n_rows; //(Modified) 20180707
}


// 2. if(fam=="poi"){
double poi_obj_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = pmin(drop(x.mat%*%b.vec),700);                                      #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>POI_BIG_XB)).fill(POI_BIG_XB);
     //   return( sum(exp(xb.vec)-y.vec*xb.vec)/length(y.vec) )
     return arma::sum(arma::exp(xb_vec)-y_vec%xb_vec)/y_vec.n_rows; //(Mod) 20180707
}

arma::vec poi_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = pmin(drop(x.mat%*%b.vec),700);                                      #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>POI_BIG_XB)).fill(POI_BIG_XB);
     //   return( drop(t(x.mat)%*%(exp(xb.vec)-y.vec))/length(y.vec) ) #
     return x_mat.t()*(arma::exp(xb_vec)-y_vec)/y_vec.n_rows;
}

arma::mat poi_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = pmin(drop(x.mat%*%b.vec),700);                                   # Mod 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>POI_BIG_XB)).fill(POI_BIG_XB);
     //   exb.vec = exp(xb.vec)                                                     # Add exb.vec
     arma::vec exb_vec = arma::exp(xb_vec); // Add exb.vec
     //   return( t(x.mat)%*%diag(exb.vec)%*%x.mat/length(y.vec))
     return x_mat.t()*arma::diagmat(exb_vec)*x_mat/y_vec.n_rows;
}

// 3. if(fam=="log"){
double log_obj_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = pmin(drop(x.mat%*%b.vec),700)                                              # Mod 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>LOGI_BIG_XB)).fill(LOGI_BIG_XB);
     //   return( sum(log(1+exp(xb.vec))-y.vec*xb.vec)/length(y.vec) )
     return arma::sum(arma::log(1+arma::exp(xb_vec))-y_vec%xb_vec)/y_vec.n_rows;
}

arma::vec log_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = pmin(drop(x.mat%*%b.vec),700)                                  #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>LOGI_BIG_XB)).fill(LOGI_BIG_XB);
     //   exb.vec = exp(xb.vec)
     arma::vec exb_vec = arma::exp(xb_vec);
     //   p.vec = exb.vec/(1+exb.vec)
     arma::vec p_vec = exb_vec/(1+exb_vec);
     //   return( drop(t(x.mat)%*%(p.vec-y.vec)/length(y.vec)) )
     return x_mat.t()*(p_vec-y_vec)/y_vec.n_rows;
}

arma::mat log_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec) {
     //   xb.vec = pmin(drop(x.mat%*%b.vec),700)                                  #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>LOGI_BIG_XB)).fill(LOGI_BIG_XB);
     //   exb.vec = exp(xb.vec)
     arma::vec exb_vec = arma::exp(xb_vec);
     //   p.vec = exb.vec/(1+exb.vec)
     arma::vec p_vec = exb_vec/(1+exb_vec);
     //   d.vec = p.vec*(1-p.vec);
     //   d.vec[d.vec<1e-7] = 1e-7
     arma::vec d_vec = p_vec%(1-p_vec);
     d_vec(arma::find(d_vec < LOGI_SMALL_NUMBER)).fill(LOGI_SMALL_NUMBER);
     //   return( t(x.mat)%*%diag(d.vec)%*%x.mat/length(y.vec) )
     return x_mat.t()*arma::diagmat(d_vec)*x_mat/y_vec.n_rows;
}

// 4. Cox
// G0.fun <- function(Y,X,S,B){
// The last element of x_mat is censored indicator.
double cox_obj_fun(arma::vec& y_vec, arma::mat& x_mat0, arma::vec& b_vec) {
     //      X <- as.matrix(X);
     // The last col of x_mat is censored indicator.
     // Rcpp::Rcout << "1" <<endl;
     arma::mat x_mat = x_mat0.cols(0, x_mat0.n_cols-2);
     // Rcpp::Rcout << "2" <<endl;
     arma::vec s_vec = x_mat0.col(x_mat0.n_cols-1);
     // Rcpp::Rcout << "3" <<endl;
     //   N <- length(Y);
     int n = y_vec.n_rows;
     //   R <- sum(S==1)
     int r = arma::sum(s_vec == 1);
     //   F.vec <- S==1
     arma::uvec f_vec = (s_vec == 1);
     //   F.mat <- matrix(rep(Y[F.vec],N),ncol=R,byrow=T) - matrix(rep(Y,R),ncol=R) <= 0
     arma::mat f_mat1(n, r);
     // Rcpp::Rcout << "4" <<endl;
     f_mat1.each_row() = y_vec(arma::find(f_vec)).t();
     // Rcpp::Rcout << "5" <<endl;

     arma::mat f_mat2(n, r);
     f_mat2.each_col() = y_vec;
     // Rcpp::Rcout << "6" <<endl;

     arma::umat f_mat = (f_mat1 - f_mat2) <= 0;
     // Rcpp::Rcout << "6.1" <<endl;

     // XB <- drop(X%*%B); XB[XB>100] <- 100
     arma::vec xb_vec = x_mat*b_vec;
     // Rcpp::Rcout << "6.2" <<endl;
     xb_vec(arma::find(xb_vec > 100)).fill(100);
     // Rcpp::Rcout << "7" <<endl;

     // Exp <- matrix(rep(exp(XB),R),ncol=R)
     arma::mat exp_mat(n, r);
     // Rcpp::Rcout << "7.1" <<endl;
     exp_mat.each_col() = arma::exp(xb_vec);
     // Rcpp::Rcout << "7.2" <<endl;

     // W0 <- pmax(colSums(F.mat*Exp),exp(-100))
     arma::vec w0_vec = arma::trans(arma::sum(f_mat%exp_mat, 0));
     // Rcpp::Rcout << "7.3" <<endl;
     w0_vec(arma::find(w0_vec < std::exp(-100))).fill(std::exp(-100));
     // Rcpp::Rcout << "8" <<endl;

     // U0 <- as.numeric(F.vec%*%XB)
     double u0 = arma::sum(f_vec%xb_vec);
     // V0 <- (-U0+sum(log(W0)))/N
     double v0 = (-u0 + arma::sum(arma::log(w0_vec)))/n;

     // return(V0)
     return v0;
}

// # score function
// G1.fun <- function(Y,X,S,B){
arma::vec cox_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat0, arma::vec& b_vec) {
     //      X <- as.matrix(X);
     // The last col of x_mat is censored indicator.
     // Rcpp::Rcout << "1" <<endl;
     arma::mat x_mat = x_mat0.cols(0, x_mat0.n_cols-2);
     // Rcpp::Rcout << "2" <<endl;
     arma::vec s_vec = x_mat0.col(x_mat0.n_cols-1);
     // Rcpp::Rcout << "3" <<endl;

     //   N <- length(Y);
     int n = y_vec.n_rows;

     //   R <- sum(S==1)
     int r = arma::sum(s_vec == 1);

     //   F.vec <- S==1
     arma::uvec f_vec = (s_vec == 1);
     //   F.mat <- matrix(rep(Y[F.vec],N),ncol=R,byrow=T) - matrix(rep(Y,R),ncol=R) <= 0
     arma::mat f_mat1(n, r);
     f_mat1.each_row() = y_vec(arma::find(f_vec)).t();

     arma::mat f_mat2(n, r);
     f_mat2.each_col() = y_vec;

     arma::umat f_mat = (f_mat1 - f_mat2) <= 0;
     //   XB <- drop(X%*%B); XB[XB>100] <- 100
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec > 100)).fill(100);

     //   Exp <- matrix(rep(exp(XB),R),ncol=R)
     arma::mat exp_mat(n, r);
     exp_mat.each_col() = arma::exp(xb_vec);

     //   U1 <- as.vector(F.vec%*%X)
     arma::vec u1_vec = x_mat.t()*f_vec;
     //   W0 <- pmax(colSums(F.mat*Exp),exp(-100))
     arma::vec w0_vec = arma::trans(arma::sum(f_mat%exp_mat, 0));
     w0_vec(arma::find(w0_vec < std::exp(-100))).fill(std::exp(-100));
     //   W1 <- t(F.mat)%*%(X*exp(XB))
     arma::mat x_exp_xb = x_mat.each_col() % arma::exp(xb_vec);
     arma::mat w1_mat = f_mat.t()*x_exp_xb;

     //      return((-U1+colSums(W1/W0))/N)
     arma::mat w1_w0 = w1_mat.each_col() / w0_vec;
     arma::vec temp = arma::trans(arma::sum(w1_w0, 0));

     return (-u1_vec + temp)/n;
}

// W1 <- t(apply(W1,1,Outer))
arma::mat apply_row_outer(arma::mat& x) {
     arma::mat ret_mat(x.n_cols*x.n_cols, x.n_rows);
     //Rcpp::Rcout << "3.1" <<endl;
     for(arma::uword i = 0; i<x.n_rows; i++) {
          //Rcpp::Rcout << "3.2" <<endl;
          arma::mat one_outer = x.row(i).t()*x.row(i);
          //Rcpp::Rcout << "3.3: one_outer.n_elem: " << one_outer.n_elem << ", ret_mat.n_rows: " << ret_mat.n_rows <<endl;
          one_outer.reshape(one_outer.n_elem, 1);
          ret_mat.col(i) = one_outer;
          //Rcpp::Rcout << "3.4" <<endl;
     }

     return ret_mat;

}
//# full hessian matrix
arma::mat cox_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat0, arma::vec& b_vec) {
     //G2.fun <- function(Y,X,S,B){
     // X <- as.matrix(X);
     // The last col of x_mat is censored indicator.
     //Rcpp::Rcout << "1" <<endl;
     arma::mat x_mat = x_mat0.cols(0, x_mat0.n_cols-2);
     arma::vec s_vec = x_mat0.col(x_mat0.n_cols-1);
     // N <- length(Y);
     int n = y_vec.n_rows;
     // R <- sum(S==1)
     int r = arma::sum(s_vec == 1);
     // F.vec <- S==1
     arma::uvec f_vec = (s_vec == 1);
     // F.mat <- matrix(rep(Y[F.vec],N),ncol=R,byrow=T) - matrix(rep(Y,R),ncol=R) <= 0
     arma::mat f_mat1(n, r);
     f_mat1.each_row() = y_vec(arma::find(f_vec)).t();

     arma::mat f_mat2(n, r);
     f_mat2.each_col() = y_vec;

     arma::umat f_mat = (f_mat1 - f_mat2) <= 0;

     //Rcpp::Rcout << "2" <<endl;
     // XB <- drop(X%*%B); XB[XB>100] <- 100
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec > 100)).fill(100);

     // Exp <- matrix(rep(exp(XB),R),ncol=R)
     arma::mat exp_mat(n, r);
     exp_mat.each_col() = arma::exp(xb_vec);

     // W0 <- pmax(colSums(F.mat*Exp),exp(-100))
     arma::vec w0_vec = arma::trans(arma::sum(f_mat%exp_mat, 0));
     w0_vec(arma::find(w0_vec < std::exp(-100))).fill(std::exp(-100));


     // W1 <- t(F.mat)%*%(X*exp(XB))
     arma::mat exp_xb = arma::exp(xb_vec);
     arma::mat x_exp_xb = x_mat.each_col() % exp_xb;
     arma::mat w1_mat = f_mat.t()*x_exp_xb;
     arma::mat w2_mat;
     //Rcpp::Rcout << "3" <<endl;
     // if(length(B)==1){
     if(b_vec.n_rows == 1) {
          //   W1 <- W1^2;
          w1_mat = arma::square(w1_mat);
          //   W2 <- t(F.mat)%*%(X^2*exp(XB))
          arma::mat x2_mat = arma::square(x_mat);
          arma::mat x2_exp_xb = x2_mat.each_col() % exp_xb;
          w2_mat = f_mat.t()*x2_exp_xb;
     } else {
          //   W1 <- t(apply(W1,1,Outer));
          // Rcpp::Rcout << "3.01" <<endl;
          w1_mat = arma::trans(apply_row_outer(w1_mat));
          //Rcpp::Rcout << "3.02" <<endl;

          //   W2 <- t(F.mat) %*% ( t( apply(X,1,Outer) ) * exp(XB) )
          arma::mat temp_mat = arma::trans(apply_row_outer(x_mat));
          w2_mat = f_mat.t()*(temp_mat.each_col()%exp_xb);
          //Rcpp::Rcout << "3.03" <<endl;
     }
     //Rcpp::Rcout << "4" <<endl;
     // return(matrix(colSums(W2/W0-W1/(W0^2)),nrow=dim(X)[2])/N)
     arma::mat ret_mat = arma::sum(w2_mat.each_col()/w0_vec - w1_mat.each_col()/arma::square(w0_vec), 0)/n;
     ret_mat.reshape(ret_mat.n_elem/x_mat.n_cols, x_mat.n_cols);
     //Rcpp::Rcout << "5" <<endl;

     return ret_mat;
}

// 5. Multinomial
// utility function for multinomial
arma::uvec y_vec_to_sy(arma::vec& y_vec) {
     arma::uword k = arma::max(y_vec);
     arma::mat y_mat(y_vec.n_rows, k-1); // remove last one as base line

     // categories
     arma::mat k_mat(1, k-1);
     for(arma::uword i = 0; i<k-1; i++) {
          k_mat(0, i) = i+1;
     }

     // for(int i = 0; i<k-1; i++) {
     //      y_mat.col(i) = y_vec;
     // }
     y_mat.each_col() = y_vec;

     arma::umat ret_mat(y_vec.n_rows, k-1);
     for(arma::uword i = 0; i<y_vec.n_rows; i++) {
          ret_mat.row(i) = (y_mat.row(i) == k_mat);
     }

     ret_mat.reshape(ret_mat.n_elem, 1);

     return ret_mat;
}

// obj.fun = function(y.vec,x.mat,b.vec){ # if y: k class = > b: p*(k-1)
double mtn_obj_fun(arma::vec& y_vec, arma::mat& sx_mat, arma::vec& sb_vec) {
     //     Rcpp::Rcout << "multinomial: v2" << endl;

     //   k = max(y.vec)
     arma::uword k = arma::max(y_vec);
     //n = length(y.vec)
     //arma::uword n = y_vec.n_rows;
     //sy.vec = as.numeric(rep(y.vec,k-1) == rep(1:(k-1),each=n))
     arma::uvec sy_vec = y_vec_to_sy(y_vec);

     //sxb.vec = pmin(drop(sx.mat%*%sb.vec),700);
     arma::vec sxb_vec = sx_mat*sb_vec;
     sxb_vec(arma::find(sxb_vec>700)).fill(700);
     //esxb.vec = exp(sxb.vec);
     arma::vec esxb_vec = arma::exp(sxb_vec);
     //esxb.mat = matrix(esxb.vec,ncol=k-1)
     arma::mat esxb_mat(esxb_vec);
     esxb_mat.reshape(esxb_mat.n_elem/(k-1), k-1);

     //loss = sum(log(1+rowSums(esxb.mat)))-sum(sy.vec*sxb.vec)
     double loss = arma::accu(arma::log(1+arma::sum(esxb_mat, 1)))-arma::accu(sy_vec%sxb_vec);
     //   return(loss/length(y.vec))
     return loss/y_vec.n_rows;
}

//obj.grad.fun = function(y.vec,x.mat,b.vec){
arma::vec mtn_obj_grad_fun(arma::vec& y_vec, arma::mat& sx_mat, arma::vec& sb_vec) {
     // k = max(y.vec)
     arma::uword k = arma::max(y_vec);
     //n = length(y.vec)
     //arma::uword n = y_vec.n_rows;
     //sy.vec = as.numeric(rep(y.vec,k-1) == rep(1:(k-1),each=n))
     arma::uvec sy_vec = y_vec_to_sy(y_vec);

     //sxb.vec = pmin(drop(sx.mat%*%sb.vec),700);
     arma::vec sxb_vec = sx_mat*sb_vec;
     sxb_vec(arma::find(sxb_vec>700)).fill(700);
     //esxb.vec = exp(sxb.vec);
     arma::vec esxb_vec = arma::exp(sxb_vec);
     //esxb.mat = matrix(esxb.vec,ncol=k-1)
     arma::mat esxb_mat(esxb_vec);
     esxb_mat.reshape(esxb_mat.n_elem/(k-1), k-1);

     //sp.mat = esxb.mat/(1+rowSums(esxb.mat));
     arma::mat sp_mat = esxb_mat.each_col()/(1+arma::sum(esxb_mat,1));
     //sp.vec = as.vector(sp.mat)
     arma::vec sp_vec = arma::vectorise(sp_mat);
     //grad.vec = -as.vector(t(sx.mat)%*%(sy.vec-sp.vec))
     arma::vec grad_vec = -arma::vectorise(sx_mat.t()*(sy_vec-sp_vec));

     //return(grad.vec/length(y.vec));
     return grad_vec/y_vec.n_rows;
}

// obj.hess.fun = function(y.vec,x.mat,b.vec){
arma::mat mtn_obj_hess_fun(arma::vec& y_vec, arma::mat& sx_mat, arma::vec& sb_vec) {
     // k = max(y.vec)
     arma::uword k = arma::max(y_vec);
     //n = length(y.vec)
     arma::uword n = y_vec.n_rows;
     //sy.vec = as.numeric(rep(y.vec,k-1) == rep(1:(k-1),each=n))
     arma::uvec sy_vec = y_vec_to_sy(y_vec);

     //sxb.vec = pmin(drop(sx.mat%*%sb.vec),700);
     arma::vec sxb_vec = sx_mat*sb_vec;
     sxb_vec(arma::find(sxb_vec>700)).fill(700);
     //esxb.vec = exp(sxb.vec);
     arma::vec esxb_vec = arma::exp(sxb_vec);
     //esxb.mat = matrix(esxb.vec,ncol=k-1)
     arma::mat esxb_mat(esxb_vec);
     esxb_mat.reshape(esxb_mat.n_elem/(k-1), k-1);
     //sp.mat = esxb.mat/(1+rowSums(esxb.mat));
     arma::mat sp_mat = esxb_mat.each_col()/(1+arma::sum(esxb_mat,1));

     //sp.vec = as.vector(sp.mat)
     arma::vec sp_vec = arma::vectorise(sp_mat);
     //sp.vec[sp.vec<1e-7] = 1e-7;
     sp_vec(arma::find(sp_vec<LOGI_SMALL_NUMBER)).fill(LOGI_SMALL_NUMBER);
     //sp.vec[sp.vec>(1-(1e-7))] = 1-(1e-7)
     sp_vec(arma::find(sp_vec>(1-LOGI_SMALL_NUMBER))).fill(1-LOGI_SMALL_NUMBER);
     //d.mat = diag(sp.vec) - outer(sp.vec,sp.vec) *  kronecker( outer(rep(1,k-1),rep(1,k-1)) , diag(n) )

     arma::mat d_mat = arma::diagmat(sp_vec) - sp_vec*sp_vec.t() %  arma::kron( arma::ones<arma::mat>(k-1, k-1) , arma::diagmat(arma::ones<arma::vec>(n)) );
     //hess.mat = t(sx.mat)%*%d.mat%*%sx.mat
     arma::mat hess_mat = sx_mat.t()*d_mat*sx_mat;

     //return(hess.mat/n)
     return hess_mat/n;
}



obj_fun_ptr get_obj_fun_ptr(std::string name) {
     if(name.compare("gaussian") == 0) {
          return lin_obj_fun;
     } else if(name.compare("poisson") == 0) {
          return poi_obj_fun;
     } else if(name.compare("binomial") == 0) {
          return log_obj_fun;
     } else if(name.compare("cox") == 0) {
          return cox_obj_fun;
     } else if(name.compare("multinomial") == 0) {
          return mtn_obj_fun;
     } else {
          throw std::invalid_argument("Invalid oject function family. Only supports \"gaussian\" (linear), \"binomial\" (logistic), \"poisson\", \"cox\", \"multinomial\".");
          return NULL;
     }
}

obj_grad_fun_ptr get_obj_grad_fun_ptr(std::string name) {
     if(name.compare("gaussian") == 0) {
          return lin_obj_grad_fun;
     } else if(name.compare("poisson") == 0) {
          return poi_obj_grad_fun;
     } else if(name.compare("binomial") == 0) {
          return log_obj_grad_fun;
     } else if(name.compare("cox") == 0) {
          return cox_obj_grad_fun;
     } else if(name.compare("multinomial") == 0) {
          return mtn_obj_grad_fun;
     } else {
          throw std::invalid_argument("Invalid oject function family @ gradient. Only supports \"gaussian\" (linear), \"binomial\" (logistic), \"poisson\", \"cox\", \"multinomial\".");
          return NULL;
     }
}

obj_hess_fun_ptr get_obj_hess_fun_ptr(std::string name) {
     if(name.compare("gaussian") == 0) {
          return lin_obj_hess_fun;
     } else if(name.compare("poisson") == 0) {
          return poi_obj_hess_fun;
     } else if(name.compare("binomial") == 0) {
          return log_obj_hess_fun;
     } else if(name.compare("cox") == 0) {
          return cox_obj_hess_fun;
     } else if(name.compare("multinomial") == 0) {
          return mtn_obj_hess_fun;
     } else {
          throw std::invalid_argument("Invalid oject function family @ hessian. Only supports \"gaussian\" (linear), \"binomial\" (logistic), \"poisson\", \"cox\", \"multinomial\".");
          return NULL;
     }
}
// End of loss functions ------------------------------------------
// ###### Newton-Raphosn for nun-penalized model : 20180707
// nr.fun = function(y.vec,x.mat,iter.max,b.eps){
arma::vec nr_fun(std::string fam, arma::vec& y_vec, arma::mat& x_mat, double iter_max, double b_eps) {
     //std::string fam = Rcpp::as<std::string>(family);

     // obj.grad.fun = get.fam.fun(fam)[[2]]
     obj_grad_fun_ptr obj_grad_fun = get_obj_grad_fun_ptr(fam);
     // obj.hess.fun = get.fam.fun(fam)[[3]]
     obj_hess_fun_ptr obj_hess_fun = get_obj_hess_fun_ptr(fam);

     //   b.vec = rep(0,dim(x.mat)[2])
     arma::vec b_vec;
     if(fam == "cox") {
          b_vec = arma::zeros<arma::vec>(x_mat.n_cols-1);
     } else {
          b_vec = arma::zeros<arma::vec>(x_mat.n_cols);
     }
     //   for(iter in 1:iter.max){
     arma::uword iter = 0;
     for(;iter <iter_max; iter++) {
          // #print(b.vec)
          //        hess.mat = obj.hess.fun(y.vec,x.mat,b.vec)
          arma::mat hess_mat = obj_hess_fun(y_vec, x_mat, b_vec);
          //        grad.vec = obj.grad.fun(y.vec,x.mat,b.vec)
          arma::vec grad_vec = obj_grad_fun(y_vec, x_mat, b_vec);
          //        nb.vec = b.vec-solve(hess.mat)%*%grad.vec
          arma::vec nb_vec = b_vec-arma::inv(hess_mat)*grad_vec;
          //        if(sum(abs(nb.vec-b.vec))<b.eps) break
          if(arma::sum(arma::abs(nb_vec-b_vec))<b_eps) break;
          //        b.vec = nb.vec
          b_vec = nb_vec;
     }
     //   return(b.vec)
     return b_vec;
}

// 3. Soft-thresholding function ------------------------------------------
// soft.fun = function(est,del){
//      return(sign(est)*(abs(est)-del)*(abs(est)>del))
// }
double soft_fun(double est, double del) {
     return sign_scalar_(est)*(std::abs(est)-del)*(std::abs(est)>del);
}
// End of soft-thresholding function ------------------------------------------

// 4. Squadratic lasso function ------------------------------------------
// it minimizes t(b.vec)%*%(q.mat)%*%b.vec + t(l.vec)%*%b.vec + lam||t(w.vec)%*%b.vec||_1


// #######################################################
// ### warm => p.eff, only name change in all functions
// ### warm => p.eff, only name change in all functions
// ### warm => p.eff, only name change in all functions
// #######################################################
double get_qlasso_fun_est(int pos, arma::mat& q_mat, arma::vec& b_vec, arma::vec& l_vec) {
     //est = -(2*sum(q.mat[-pos,pos]*b.vec[-pos])+l.vec[pos])/q.mat[pos,pos]/2;

     // arma::vec q_mat_col = q_mat.col(pos);
     // return -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;

     arma::vec q_vec1 = q_mat.col(pos); q_vec1.shed_row(pos);
     arma::vec b_vec1 = b_vec; b_vec1.shed_row(pos);

     return -(2*arma::sum(q_vec1%b_vec1)+l_vec(pos))/q_mat(pos, pos)/2;
}

double get_qlasso_fun_del(int pos, arma::mat& q_mat, arma::vec& w_vec, double lam) {
     //del = lam*w.vec[pos]/q.mat[pos,pos]/2
     return lam*w_vec(pos)/q_mat(pos, pos)/2;
}

//qlasso.fun = function(q.mat,l.vec,b.vec,w.vec,lam,iter.max,b.eps,k.eps,p.eff,q.rank){
//const double SMALL_NUMBER = 1e-7;

// // [[Rcpp::export]]
// Rcpp::List qlasso_fun(arma::mat& q_mat, arma::vec& l_vec, arma::vec& b_vec0, arma::vec& w_vec,
//                                   double lam, double iter_max, double b_eps, double k_eps, arma::uword p_eff, arma::uword q_rank) {
int qlasso_fun(arma::mat& q_mat, arma::vec& l_vec, arma::vec& b_vec0, arma::vec& w_vec,
               double lam, double iter_max, double iiter_max, double b_eps, double k_eps,
               arma::uword p_eff, arma::uword q_rank, bool cut, double c_eps, p_ncpen_ret& ret_buff) {

     arma::vec b_vec = b_vec0;
     //   p = length(b.vec)
     //   arma::uword p = b_vec.n_rows;
     //   #ob.vec = b.vec                                   ### (Del)

     //Rcout << "qlasso_fun v1.\n";
     //f.vec = rep(0,iter.max)
     arma::vec f_vec = arma::zeros<arma::vec>(iter_max);

     // Vriables for loop resue----------------------
     arma::vec ob_vec = b_vec;
     //arma::vec g_vec;// = arma::zeros<arma::vec>(p);
     //arma::vec p_vec;// = arma::zeros<arma::vec>(p);
     arma::vec g_vec(b_vec.n_rows);
     arma::vec p_vec(b_vec.n_rows);
     bool kkt0 = false;
     bool kkt1 = false;
     // Variables for loop resue----------------------
     arma::uword iter = 0;
     arma::uvec a_set;
     //   for(iter in 1:iter.max){#iter
     for(;iter < iter_max; iter++){ // iter
          // ##### active set
          //p0.eff = p.eff                                  ### (Add)
          arma::uword pp_eff = p_eff; //### p0.eff => pp.eff

          //a1.set = c(1:p)[(b.vec!=0)&(w.vec!=0)]
          arma::uvec a1_set = arma::find( ((b_vec!= 0) % (w_vec!=0))> 0 );
          arma::uvec a2_set = arma::find(w_vec == 0);
          a_set = arma::join_cols(a1_set, a2_set);
          //for(iiter in 1:iter.max){#iiter
          //for(arma::uword iiter = 0; iiter < iter_max; iiter++) { // iiter
          arma::uword iiter = 0;
          bool con = false;
          for(; iiter < iiter_max; iiter++) { // iiter
               ob_vec = b_vec;
               //ab.vec = b.vec[a.set]
               arma::vec ab_vec = b_vec(a_set);
               //for(pos in a1.set){
               for(arma::uword idx = 0; idx < a1_set.n_rows; idx++) {
                    arma::uword pos = a1_set(idx);
                    //est = -(2*sum(q.mat[-pos,pos]*b.vec[-pos])+l.vec[pos])/q.mat[pos,pos]/2
                    //double est = -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;
                    double est = get_qlasso_fun_est(pos, q_mat, b_vec, l_vec);

                    //del = lam*w.vec[pos]/q.mat[pos,pos]/2
                    //double del = lam*w_vec(pos)/q_mat(pos, pos)/2;
                    double del = get_qlasso_fun_del(pos, q_mat, w_vec, lam);

                    //b.vec[pos] = soft.fun(est,del)
                    b_vec(pos) = soft_fun(est, del);
               }

               //for(pos in a2.set){
               for(arma::uword idx = 0; idx < a2_set.n_rows; idx++) {
                    arma::uword pos = a2_set(idx);

                    //est = -(2*sum(q.mat[-pos,pos]*b.vec[-pos])+l.vec[pos])/q.mat[pos,pos]/2
                    //double est = -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;
                    double est = get_qlasso_fun_est(pos, q_mat, b_vec, l_vec);

                    //b.vec[pos] = est
                    b_vec(pos) = est;
               }

               con = arma::sum(arma::abs(b_vec-ob_vec)) < b_eps;
               // if(con) break                                                   ### add con
               if(con) break;

               //#####################################
               //##### projection: reduced model #####
               //#####################################
               //b.vec[abs(b.vec)<c.eps] = 0                                ### add this
               if(cut == true) {
                    b_vec(arma::find(arma::abs(b_vec) < c_eps)).fill(0);
               }

               //############################################# Add
               //##### iiter projection
               //if(iiter>p0.eff){
               if(iiter > pp_eff) { //### p0.eff => pp.eff
                    arma::vec oob_vec = b_vec;
                    //a.set = c(1:p)[b.vec!=0]
                    a_set = arma::find(b_vec != 0);

                    //if(length(a.set)<q.rank){
                    if(a_set.n_rows <= q_rank){ //### please check the inequality here! if length(a.set)>q.rank=n then the inverse may not exist!!!
                         //pb.vec = b.vec*0
                         arma::vec pb_vec = arma::zeros<arma::vec>(b_vec.n_rows);

                         //pb_vec(a_set) = -arma::pinv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                         pb_vec(a_set) = -arma::inv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                         //pbf = t(pb.vec[a.set])%*%(q.mat[a.set,a.set])%*%pb.vec[a.set] + sum(l.vec[a.set]*pb.vec[a.set])+lam*sum(abs(w.vec[a.set]*pb.vec[a.set]))
                         double pbf = arma::as_scalar(pb_vec(a_set).t()*q_mat(a_set, a_set)*pb_vec(a_set)) + arma::sum(l_vec(a_set)%pb_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%pb_vec(a_set)));
                         //bf  = t( b.vec[a.set])%*%(q.mat[a.set,a.set])%*% b.vec[a.set] + sum(l.vec[a.set]* b.vec[a.set])+lam*sum(abs(w.vec[a.set]* b.vec[a.set]))
                         double bf = arma::as_scalar(b_vec(a_set).t()*q_mat(a_set, a_set)*b_vec(a_set)) + arma::sum(l_vec(a_set)%b_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%b_vec(a_set)));
                         //[e 2017011]

                         //if(pbf<bf){ ob.vec = b.vec = pb.vec; cat("iiter projenction works","\n")
                         if(pbf < bf){ // How about if((pbf-bf)<b_eps) {
                              //ob.vec = b.vec = pb.vec;
                              // ob_vec = pb_vec; # deleted 20180707
                              b_vec = pb_vec;
                              //cat("projection on active set in quad.lasso works.","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                              if(NCPEN_DEVELOP_MODE == true) {
                                 Rcpp::Rcout << "Projection on active set in qlasso_fun works. lam=" << lam << ", df=" << sum(b_vec!=0) << ", iiter=" << iiter << std::endl;
                              }
                              //} else { p0.eff = 2*p0.eff; cat("iiter projenction fails","\n");  }
                         } else {
                              //p0.eff = 2*p0.eff;
                              //cat("projection on active set in quad.lasso fails.","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                              if(NCPEN_DEVELOP_MODE == true) {
                                 Rcpp::Rcout << "projection on active set in qlasso_fun fails. lam" << lam << ", df=" << sum(b_vec!=0) << ", iiter=" << iiter << std::endl;
                              }

                              b_vec = oob_vec;
                         }
                         pp_eff = 2*pp_eff;
                    }
               }

               con = arma::sum(arma::abs(b_vec-ob_vec)) < b_eps;
               // if(con) break                                                   ### add con
               if(con) break;

          }//iiter

          if(NCPEN_DEVELOP_MODE == true) {
               Rcpp::Rcout << "qlasso_fun : a_set=" << a_set.n_rows << ", lambda=" << lam << ", iiter=" << iiter << ", con=" << con << std::endl;
          }

          //b.vec[abs(b.vec)<1e-7] = 0      ### (Add)
          //[s/e 20170111] b_vec(arma::find(arma::abs(b_vec) < SMALL_NUMBER)).fill(0);

          //###################
          //##### null set
          //###################
          //n1.set = c(1:p)[(b.vec==0)&(w.vec!=0)]
          arma::uvec n_set = arma::find( ((b_vec==0) % (w_vec!=0))> 0 );


          //for(pos in n1.set){
          for(arma::uword idx = 0; idx < n_set.n_rows; idx++) {
               arma::uword pos = n_set(idx);
               //est = -(2*sum(q.mat[-pos,pos]*b.vec[-pos])+l.vec[pos])/q.mat[pos,pos]/2
               //double est = -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;
               double est = get_qlasso_fun_est(pos, q_mat, b_vec, l_vec);
               //del = lam*w.vec[pos]/q.mat[pos,pos]/2
               //double del = lam*w_vec(pos)/q_mat(pos, pos)/2;
               double del = get_qlasso_fun_del(pos, q_mat, w_vec, lam);
               //b.vec[pos] = soft.fun(est,del)
               b_vec(pos) = soft_fun(est, del);
          }


          //b.vec[abs(b.vec)<1e-7] = 0    ### (Add)
          if(cut == true) {
               b_vec(arma::find(arma::abs(b_vec) < c_eps)).fill(0);
          }


          // ##################################
          // ##### projection: full model #####
          // ##################################
          //if(iter>p.eff){
          if(iter > p_eff){ //### please check the inequality here! if length(a.set)>q.rank=n then the inverse may not exist!!!
               arma::vec oob_vec = b_vec;
               //a.set = c(1:p)[b.vec!=0]
               arma::uvec a_set = arma::find(b_vec != 0);
               //if(length(a.set)<=q.rank){
               if(a_set.n_rows <= q_rank) {
                    //pb.vec = b.vec*0
                    arma::vec pb_vec = b_vec*0;

                    //pb_vec(a_set) = -arma::pinv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                    pb_vec(a_set) = -arma::inv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                    ///pbf = t(pb.vec[a.set])%*%(q.mat[a.set,a.set])%*%pb.vec[a.set] + sum(l.vec[a.set]*pb.vec[a.set])+lam*sum(abs(w.vec[a.set]*pb.vec[a.set]))
                    double pbf = arma::as_scalar(pb_vec(a_set).t()*q_mat(a_set, a_set)*pb_vec(a_set)) + arma::sum(l_vec(a_set)%pb_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%pb_vec(a_set)));
                    //bf  = t( b.vec[a.set])%*%(q.mat[a.set,a.set])%*% b.vec[a.set] + sum(l.vec[a.set]* b.vec[a.set])+lam*sum(abs(w.vec[a.set]* b.vec[a.set]))
                    double bf = arma::as_scalar(b_vec(a_set).t()*q_mat(a_set, a_set)*b_vec(a_set)) + arma::sum(l_vec(a_set)%b_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%b_vec(a_set)));
                    //[e 20170111]

                    if(pbf<bf){
                         //b.vec = pb.vec;
                         b_vec = pb_vec;
                         //cat("projection on full set in quad.lasso works.","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                         if(NCPEN_DEVELOP_MODE == true) {
                              Rcpp::Rcout << "Projection on full set in qlasso_fun works. lam=" << lam << ", df=" << sum(b_vec!=0) << ", iter=" << iter << std::endl;
                         }
                    } else {
                         b_vec = oob_vec;
                         //cat("projection on full set in quad.lasso fails.","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                         if(NCPEN_DEVELOP_MODE == true) {
                              Rcpp::Rcout << "Projection on full set in qlasso_fun fails. lam" << lam << ", df=" << sum(b_vec!=0) << ", iter=" << iter << std::endl;
                         }
                         //Rcout << "projection on full set in quad.lasso fails. lam" << lam << ", df=" << sum(b_vec!=0) << std::endl;
                    }
                    //p.eff = 2*p.eff;
                    p_eff = 2*p_eff;
               }
          }
          //##### iter projection

          //b.vec[abs(b.vec)<b.eps] = 0                                ### add this
          if(cut == true) {
               b_vec(arma::find(arma::abs(b_vec) < c_eps)).fill(0);
          }


          // #####################
          // ##### check KKT #####
          // #####################
          //f.vec[iter] = t(b.vec)%*%q.mat%*%b.vec+sum(l.vec*b.vec)+lam*sum(abs(w.vec*b.vec))
          f_vec(iter) = arma::as_scalar(b_vec.t()*q_mat*b_vec)+arma::sum(l_vec%b_vec) +lam*arma::sum(arma::abs(w_vec%b_vec));

          //Moved to inside the for loop, was outside
          //g.vec = 2*drop(q.mat%*%b.vec)+l.vec
          g_vec = 2*q_mat*b_vec+l_vec;
          //p.vec = lam*w.vec*sign(b.vec)
          p_vec = lam*w_vec%arma::sign(b_vec);
          a_set = arma::find( ((b_vec!=0) + (w_vec==0)) >0 ); //                                               #### (Add
          n_set = arma::find( ((b_vec==0) % (w_vec!=0)) >0 ); //                                               #### (Add
          //kkt0 = sum(abs((g.vec[a.set]+p.vec[a.set]))<k.eps)==length(a.set)      #### (Mod) sum => sum(abs( )); last sum => length
          kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps)==a_set.n_rows;    //  #### (Mod) sum => sum(abs( )); last sum => length
          //kkt0 = get_kkt0(b_vec, g_vec, p_vec, k_eps);
          //kkt1 = sum(abs(g.vec[n.set])-lam*w.vec[n.set]<k.eps)==length(n.set)    #### (Mod) Check parentheses ¸Last· sum => length
          kkt1 = arma::sum(arma::abs(g_vec(n_set))-lam*w_vec(n_set)<k_eps)==n_set.n_rows;
          //kkt1 = get_kkt1(b_vec, w_vec, g_vec, p_vec, lam, k_eps);
          //Moved to inside the for loop, was outside

          //#if(sum(abs(b.vec-ob.vec))<b.eps) break    ### (Del)
          //if(kkt0&kkt1) break #iter break            ### (Add)
          if(kkt0&&kkt1) {
               iter++; // for f_vec.head(iter);
               break; // iter break                                    ### (Add)
          }
          // #ob.vec = b.vec                            ### (Del)
     }//#iter

     if(NCPEN_DEVELOP_MODE == true) {
          // bool con2 = kkt0&&kkt1;
          Rcpp::Rcout << "qlasso_fun : a_set=" << a_set.n_rows << ", lambda=" << lam << ", iter=" << iter << ", kkt0=" << kkt0 << ", kkt1=" << kkt1 << std::endl;
     }
     //return(list(g.vec=g.vec+p.vec,b.vec=b.vec,f.vec=f.vec[1:iter],con=(kkt0&kkt1)))
     arma::vec f_vec_head = f_vec.head(iter);
     // return Rcpp::List::create(Rcpp::Named("g.vec") = g_vec + p_vec,
     //                           Rcpp::Named("b.vec") = b_vec,
     //                           Rcpp::Named("f.vec") = f_vec_head,
     //                           Rcpp::Named("con") = kkt0 && kkt1);
     ret_buff.g_vec = g_vec + p_vec;
     ret_buff.b_vec = b_vec;
     ret_buff.f_vec = f_vec_head;
     ret_buff.con = kkt0 && kkt1;

     // std::cout << "New NCPEN is Running...." << endl;

     return 0;
}

// pointwise nonconvex penalized estimation
// p.ncpen.fun = function(y.vec,x.mat,b.vec,w.vec,lam,gam,tau,iter.max,b.eps,k.eps,p.eff){
int p_ncpen_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec0, arma::vec& w_vec,
                double lam, double gam, double tau, double alp, double iter_max, double qiter_max, double qiiter_max,
                double b_eps, double k_eps, arma::uword p_eff, bool cut, double c_eps,
                obj_fun_ptr obj_fun, obj_grad_fun_ptr obj_grad_fun, obj_hess_fun_ptr obj_hess_fun,
                pen_fun_ptr pen_fun, pen_grad_fun_ptr pen_grad_fun, p_ncpen_ret& ret_buff) {

     //q.rank = dim(x.mat)[2]        ### (Mod)
     int q_rank = x_mat.n_cols;
     //p = length(b.vec)             ### (Add)
     arma::vec b_vec = b_vec0;
     //int p = b_vec.n_rows;
     //ob.vec = b.vec
     arma::vec ob_vec = b_vec;

     arma::vec g_vec(b_vec.n_rows);
     arma::vec p_vec(b_vec.n_rows);
     arma::mat q_mat(b_vec.n_rows, b_vec.n_rows);

     bool kkt0 = false;
     bool kkt1 = false;
     //f.vec = rep(0,iter.max)
     arma::vec f_vec = zeros<vec>(iter_max);

     arma::uvec a_set;

     //for(iter in 1:iter.max){
     arma::uword iter = 0;
     for(; iter<iter_max; iter++) {
          //h.mat = obj.hess.fun(y.vec,x.mat,b.vec)
          arma::mat h_mat = obj_hess_fun(y_vec, x_mat, b_vec);
          //g.vec = obj.grad.fun(y.vec,x.mat,b.vec)
          g_vec = obj_grad_fun(y_vec, x_mat, b_vec);
          //q.mat = h.mat/2 + (1-alp)*lam*w.vec*diag(rep(1,p))
          arma::mat w_diag_mat = arma::zeros<mat>(b_vec.n_rows, b_vec.n_rows);
          w_diag_mat.diag() = w_vec;
          q_mat = h_mat/2 + (1-alp)*lam*w_diag_mat;

          //l.vec = g.vec - drop(h.mat%*%b.vec) + alp*w.vec*pen.grad.fun(b.vec,lam,gam,tau) - alp*lam*w.vec*sign(b.vec)
          arma::vec l_vec = g_vec - h_mat*b_vec + alp*w_vec%pen_grad_fun(b_vec, lam, gam, tau) - alp*lam*w_vec%arma::sign(b_vec);

          //if(pen=="sridge"){ q.mat = q.mat + gam*w.vec*diag(rep(1,p))/2; l.vec = l.vec - gam*w.vec*b.vec }
          if(pen_fun == sridge_pen_fun) {
               q_mat = q_mat + gam*w_diag_mat/2;
               l_vec = l_vec - gam*w_vec%b_vec;
               // Rcpp::Rcout << "pen_fun == sridge_pen_fun works." << endl;
          }

          //b.vec = qlasso.fun(q.mat,l.vec,b.vec,w.vec,lam*alp,iter.max,b.eps,k.eps,p.eff,q.rank)$b.vec
          p_ncpen_ret ret_buff;
          qlasso_fun(q_mat, l_vec, b_vec, w_vec, lam*alp, qiter_max, qiiter_max, b_eps, k_eps, p_eff, q_rank, cut, c_eps, ret_buff);
          b_vec = ret_buff.b_vec;

          //t.vec = w.vec*pen.grad.fun(ob.vec,lam,gam,tau)-lam*w.vec*sign(ob.vec)                          ### (Add)  use ob.vec
          arma::vec t_vec = w_vec%(pen_grad_fun(ob_vec, lam, gam, tau)-lam*arma::sign(ob_vec)); // (Add)  use ob.vec

          //if(pen=="sridge"){ t.vec = t.vec - gam*w.vec*ob.vec }
          if(pen_fun == sridge_pen_fun){ t_vec = t_vec - gam*w_vec%ob_vec; }

          //oob = obj.fun(y.vec,x.mat,ob.vec,r.eff) + sum(t.vec*ob.vec) + lam*sum(w.vec*abs(ob.vec))   ### (Add)  use ob.vec
          //double oob = obj_fun(y_vec, x_mat, ob_vec, r_eff) + arma::sum(t_vec%ob_vec) +     lam*arma::sum(w_vec%arma::abs(ob_vec)); // (Add)  use ob.vec
          //oob = obj.fun(y.vec,x.mat,ob.vec) + alp*sum(t.vec*ob.vec) + alp*lam*sum(w.vec*abs(ob.vec)) + (1-alp)*lam*sum(w.vec*ob.vec^2)
          double oob = obj_fun(y_vec, x_mat, ob_vec) + alp*arma::sum(t_vec%ob_vec) + alp*lam*arma::sum(w_vec%arma::abs(ob_vec)) + (1-alp)*lam*arma::sum(w_vec%arma::pow(ob_vec, 2));

          //ob  = obj.fun(y.vec,x.mat, b.vec,r.eff) + sum(t.vec* b.vec) + lam*sum(w.vec*abs( b.vec))   ### (Add)  use b.vec
          //double ob  = obj_fun(y_vec, x_mat,  b_vec, r_eff) + arma::sum(t_vec% b_vec) + lam*arma::sum(w_vec%arma::abs( b_vec)); // (Add)  use b.vec
          //ob  = obj.fun(y.vec,x.mat, b.vec) + alp*sum(t.vec* b.vec) + alp*lam*sum(w.vec*abs( b.vec)) + (1-alp)*lam*sum(w.vec* b.vec^2)
          double ob = obj_fun(y_vec, x_mat, b_vec) + alp*arma::sum(t_vec%b_vec) + alp*lam*arma::sum(w_vec%arma::abs(b_vec)) + (1-alp)*lam*arma::sum(w_vec%arma::pow(b_vec, 2));

          //if(pen=="sridge"){ oob = oob + gam*sum(w.vec*ob.vec^2)/2; ob = ob + gam*sum(w.vec*b.vec^2)/2 }
          if(pen_fun == sridge_pen_fun) {
               oob = oob + gam*arma::sum(w_vec%arma::pow(ob_vec,2))/2;
               ob  =  ob + gam*arma::sum(w_vec%arma::pow(b_vec,2))/2;
          }


          //#################################################### if - for Add
          //#################################################### if - for Add
          //if(ob>oob){#mlqa
          //if(ob>oob) { // mlqa
          //if(ob>(oob+b.eps)){#mlqa
          if(ob>oob+b_eps) {
               //cat("starting golden section algorithm in pncpen.fun","\n")
               if(NCPEN_DEVELOP_MODE == true) {
                    Rcpp::Rcout << "starting golden section algorithm in p_ncpen_fun." << std::endl;
               }
               //for(iiter in 1:iter.max){
               //for(arma::uword iiter = 0; iiter < iter_max; iiter++) { // iiter

               //cat("start golden section algorithm in pncpen.fun","\n")
               //gold = (sqrt(5)-1)/2;a = -1; b = 1
               double gold = (std::sqrt(5)-1)/2; double a = -1; double b = 1;
               arma::uword iiter = 0;
               for(; iiter < iter_max; iiter++) { // iiter
                    //c = gold*a+(1-gold)*b;
                    double c = gold*a+(1-gold)*b;
                    //d = (1-gold)*a+gold*b
                    double d = (1-gold)*a+gold*b;


                    //b1.vec = c*b.vec + (1-c)*ob.vec;
                    arma::vec b1_vec = c* b_vec + (1-c)*ob_vec;
                    //b2.vec = d*b.vec + (1-d)*ob.vec
                    arma::vec b2_vec = d*b_vec + (1-d)*ob_vec;

                    //ob1 = obj.fun(y.vec,x.mat,b1.vec) + alp*sum(t.vec*b1.vec) + alp*lam*sum(w.vec*abs(b1.vec)) + (1-alp)*lam*sum(w.vec*b1.vec^2)
                    double ob1 = obj_fun(y_vec, x_mat, b1_vec) + alp*arma::sum(t_vec%b1_vec) + alp*lam*arma::sum(w_vec%arma::abs(b1_vec)) + (1-alp)*lam*arma::sum(w_vec%arma::pow(b1_vec, 2));

                    //ob2 = obj.fun(y.vec,x.mat,b2.vec) + alp*sum(t.vec*b2.vec) + alp*lam*sum(w.vec*abs(b2.vec)) + (1-alp)*lam*sum(w.vec*b2.vec^2)
                    double ob2 = obj_fun(y_vec, x_mat, b2_vec) + alp*arma::sum(t_vec%b2_vec) + alp*lam*arma::sum(w_vec%arma::abs(b2_vec)) + (1-alp)*lam*arma::sum(w_vec%arma::pow(b2_vec, 2));

                    //if(pen=="sridge"){ ob1 = ob1 + gam*sum(w.vec*b1.vec^2)/2; ob2 = ob2 + gam*sum(w.vec*b2.vec^2)/2 }
                    if(pen_fun == sridge_pen_fun) {
                         ob1 = ob1 + gam*arma::sum(w_vec%arma::pow(b1_vec,2))/2;
                         ob2 = ob2 + gam*arma::sum(w_vec%arma::pow(b2_vec,2))/2;
                    }

                    //#print(c(ob1,ob2,abs(ob1-ob2)))
                    //if(ob1 > ob2){ a = c; } else { b = d; }   #### Mod
                    if(ob1 > ob2){ a = c; } else { b = d; }   // Mod
                    //if(sum(abs(b.vec-ob.vec))<b.eps) break
                    if(arma::sum(arma::abs(b_vec-ob_vec))<b_eps) { break; }
               }
               //Rcout << "CPP p_ncpen_fun | iiter:" << iiter << std::endl;
          } //mlqa
          //#################################################### if - for Add
          //#################################################### if - for Add

          //f.vec[iter] = obj.fun(y.vec,x.mat,b.vec) + alp*sum(w.vec*pen.fun(b.vec,lam,gam,tau)) + (1-alp)*lam*sum(w.vec*b.vec^2)
          f_vec(iter) = obj_fun(y_vec, x_mat, b_vec) + alp*arma::sum(w_vec%pen_fun(b_vec, lam, gam, tau)) + (1-alp)*lam*arma::sum(w_vec%arma::pow(b_vec,2));


          //########## moved to inside for loop; was outside
          //g.vec = obj.grad.fun(y.vec,x.mat,b.vec) + alp*w.vec*pen.grad.fun(b.vec,lam,gam,tau) - alp*lam*w.vec*sign(b.vec) + 2*(1-alp)*lam*w.vec*b.vec
          g_vec = obj_grad_fun(y_vec, x_mat, b_vec) + alp*w_vec%pen_grad_fun(b_vec, lam, gam, tau) - alp*lam*w_vec%arma::sign(b_vec) + 2*(1-alp)*lam*w_vec%b_vec;



          //p.vec = alp*lam*w.vec*sign(b.vec);
          p_vec = alp*lam*w_vec%arma::sign(b_vec);

          a_set = arma::find( ((b_vec!=0) + (w_vec==0)) >0 ); //                                               #### Add
          arma::uvec n_set = arma::find( ((b_vec==0) % (w_vec!=0)) >0 ); //                                               #### (Add

          //kkt0 = sum(abs(g.vec[a.set]+p.vec[a.set])<k.eps)==length(a.set)     # Mod sum(abs( )<k.eps)==  a.set; last sum => length
          kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps) == a_set.n_rows;
          //kkt1 = sum(abs(g.vec[n.set])-alp*lam*w.vec[n.set]<k.eps)==length(n.set)
          kkt1 = arma::sum(arma::abs(g_vec(n_set))-alp*lam*w_vec(n_set)<k_eps)== n_set.n_rows;
          //kkt0 = get_kkt0(b_vec, g_vec, p_vec, k_eps);
          //kkt1 = get_kkt1(b_vec, w_vec, g_vec, p_vec, lam, k_eps);
          //########## moved to inside for loop; was outside

          //#if(sum(abs(b.vec-ob.vec))<b.eps) break                            ### (Del)
          if(NCPEN_DEVELOP_MODE == true) {
               //bool con2 = kkt0&&kkt1;
               Rcpp::Rcout << "_________p_ncpen_fun : a_set=" << a_set.n_rows <<
                    ", lambda=" << lam << ", iter=" << iter << ", kkt0=" << kkt0 << ", kkt1=" << kkt1 << std::endl;
          }

          //if(kkt0&kkt1) break #iter break                                    ### (Add)
          if(kkt0&&kkt1) {
               iter++; // for f_vec.head(iter);
               break; // iter break                                    ### (Add)
          }
          //ob.vec = b.vec
          ob_vec = b_vec;
     }

     if(NCPEN_DEVELOP_MODE == true) {
          //bool con2 = kkt0&&kkt1;
          Rcpp::Rcout << "_________p_ncpen_fun : a_set=" << a_set.n_rows <<
               ", lambda=" << lam << ", iter=" << iter << ", kkt0=" << kkt0 << ", kkt1=" << kkt1 << std::endl;
     }

     //return(list(g.vec=g.vec+p.vec,b.vec=b.vec,f.vec=f.vec[1:iter],con=(kkt0&kkt1)))
     arma::vec f_vec_head = f_vec.head(iter);
     // return Rcpp::List::create(Rcpp::Named("g.vec") = g_vec + p_vec,
     //                           Rcpp::Named("b.vec") = b_vec,
     //                           Rcpp::Named("f.vec") = f_vec_head,
     //                           Rcpp::Named("con") = kkt0 && kkt1);

     ret_buff.g_vec = g_vec + p_vec;
     ret_buff.b_vec = b_vec;
     ret_buff.f_vec = f_vec_head;
     ret_buff.con = kkt0 && kkt1;

     return 0;
}
// end of pointwise nonconvex penalized estimation

//------------------------------------------------
// nonconvex penalized estimation
// r.eff Add
// p.max moved
// warm to p.eff

// ncpen.fun = function(y.vec,x.mat,w.vec,lam.vec,gam,tau,alp,
//                      d.max,iter.max,b.eps,k.eps,p.eff,fam,pen,loc,ob.vec,div){

int ncpen_fun(arma::vec& y_vec, arma::mat& x_mat0,arma::vec& w_vec0, arma::vec& lam_vec0,
              double gam, double tau, double alp, arma::uword d_max, double iter_max, double qiter_max, double qiiter_max, double b_eps, double k_eps,
              arma::uword p_eff, bool cut, double c_eps, arma::uword add,
              std::string fam, std::string pen,
              bool loc, arma::vec& ob_vec, int div,
              ncpen_ret& ret_buff) {

     arma::uvec warning_code_vec = arma::zeros<arma::uvec>(10);

     obj_fun_ptr      obj_fun = get_obj_fun_ptr(fam);
     obj_grad_fun_ptr obj_grad_fun = get_obj_grad_fun_ptr(fam);
     obj_hess_fun_ptr obj_hess_fun = get_obj_hess_fun_ptr(fam);
     pen_fun_ptr pen_fun = get_pen_fun_ptr(pen);
     pen_grad_fun_ptr pen_grad_fun = get_pen_grad_fun_ptr(pen);

     // Create local copy of x_mat0
     arma::mat x_mat = x_mat0;
     arma::vec w_vec = w_vec0;
     arma::vec lam_vec = lam_vec0;



     // ##### mod
     // if(fam!="cox"){ p = dim(x.mat)[2] } else { p = dim(x.mat)[2]-1 }
     arma::uword p = x_mat.n_cols;
     if(fam == "cox") {p = p-1; } // last column is ceonsored indicator vector

     arma::vec b_vec(p);
     arma::uvec a_set;

     if(pen == "ridge") {
          b_vec.fill(0.001);
          a_set = linspace<arma::uvec>(0, p-1, p);
          d_max = p;
          cut = false;
     } else {
          b_vec.fill(0);
          a_set = arma::find( ((b_vec!=0) + (w_vec==0)) >0 );
          arma::mat ax_mat;
          if(a_set.n_rows > 0) {
               if(fam != "cox") {
                    ax_mat = x_mat.cols(a_set);
               } else {
                    arma::uvec ca_set = a_set;
                    ca_set.resize(ca_set.n_rows+1);
                    // index starts from 0.
                    // At the beginning, p = p-1.
                    // So, p is the index for the last columin of x_mat which contains censored indicators.
                    ca_set(ca_set.n_rows-1) = p;
                    ax_mat = x_mat.cols(ca_set);
               }
               b_vec(a_set) = nr_fun(fam, y_vec, ax_mat, iter_max, b_eps);
          }

     }

     arma::uvec n_set = arma::find( ((b_vec==0) % (w_vec!=0)) >0 );
     //r = length(lam.vec)
     arma::uword r = lam_vec.n_rows;

     //b.mat = matrix(0,p,r)
     arma::mat b_mat = arma::zeros<arma::mat>(p, r);
     //g.mat = matrix(0,p,r)
     arma::mat g_mat = arma::zeros<arma::mat>(p, r);
     //c.mat = matrix(0,2,r)
     arma::mat c_mat = arma::zeros<arma::mat>(2, r);
     //f.vec = rep(0,r)
     arma::vec f_vec = arma::zeros<arma::vec>(r);
     //d.vec = rep(0,r)
     arma::vec d_vec = arma::zeros<arma::vec>(r);

     //arma::vec g_vec(b_vec.n_rows);
     //arma::vec p_vec(b_vec.n_rows);
     bool kkt0 = false;
     bool kkt1 = false;

     //for(pos in 1:r){#pos
     arma::uword pos = 0;
     for(;pos < r; pos++) { //#pos
          //Rcout <<"for pos:" <<pos<<endl;
          //lam = lam.vec[pos]
          double lam = lam_vec(pos);

          //if(loc==TRUE) b.vec = ob.vec
          if(loc == true) {b_vec = ob_vec;}

          //g.vec = obj.grad.fun(y.vec,x.mat,b.vec) + alp*w.vec*pen.grad.fun(b.vec,lam,gam,tau) - alp*lam*w.vec*sign(b.vec) + 2*(1-alp)*lam*w.vec*b.vec
          arma::vec g_vec = obj_grad_fun(y_vec, x_mat, b_vec) + alp*w_vec%pen_grad_fun(b_vec, lam, gam, tau) - alp*lam*w_vec%arma::sign(b_vec) + 2*(1-alp)*lam*w_vec%b_vec;
          //p.vec = alp*lam*w.vec*sign(b.vec);
          arma::vec p_vec = alp*lam*w_vec%arma::sign(b_vec);

          //a.set = c(1:p)[b.vec!=0]
          a_set = arma::find(b_vec!=0);
          //n.set = c(1:p)[b.vec==0]
          arma::uvec n_set = arma::find(b_vec==0);
          //kkt0 = sum(abs(g.vec[a.set]+p.vec[a.set])<k.eps)==length(a.set)     # Mod sum(abs( )<k.eps)==  a.set; last sum => length
          kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps) == a_set.n_rows;
          //kkt1 = sum(abs(g.vec[n.set])-alp*lam*w.vec[n.set]<k.eps)==length(n.set)
          kkt1 = arma::sum(arma::abs(g_vec(n_set))-alp*lam*w_vec(n_set)<k_eps)== n_set.n_rows;

          //if(!(kkt0&kkt1)){
          arma::uword iter = 0;
          if(!(kkt0&&kkt1)){
               //for(iter in 1:iter.max){#iter
               for(; iter < iter_max; iter++) { //#iter
                    //p.set = n.set[(abs(g.vec[n.set])-alp*lam*w.vec[n.set])>k.eps]
                    arma::uvec p_set = n_set( arma::find( (arma::abs(g_vec(n_set))-alp*lam*w_vec(n_set))>k_eps ) );

                    //if(length(p.set)<=1){
                    if(p_set.n_rows == 1) {
                         //ad = p.set[which.max(abs(g.vec[p.set]))]
                         arma::vec abs_g_p = arma::abs(g_vec(p_set));
                         if(abs_g_p.n_rows > 0) {
                              //arma::uword ad = p_set(which_max(abs_g_p)); // ad is scalar
                              arma::uword ad = p_set(arma::index_max(abs_g_p)); // ad is scalar
                              //a.set = c(a.set,ad) ### mod add => ad
                              a_set.resize(a_set.n_rows + 1);
                              a_set(a_set.n_rows-1) = ad;
                         }
                    } else if (p_set.n_rows > 1) {
                         //od = order(abs(g.vec[p.set]),decreasing=T)
                         arma::vec abs_g_p = arma::abs(g_vec(p_set));
                         arma::uvec od = sort_index(abs_g_p, "descend");
                         //ed = min(length(p.set),add)
                         arma::uword ed = std::min(p_set.n_rows, add);
                         //ad= p.set[od[1:ed]]
                         arma::uvec ad = p_set(od.rows(0, ed-1));
                         //a.set = c(a.set,ad) ### mod add => ad
                         a_set = join_cols(a_set, ad); //### mod add => ad
                    }

                    //#############################################################################################################
                    //### R reconize 1 column matrix as vector.
                    // ax.mat = matrix(x.mat[,a.set],ncol=length(a.set))
                    // b.vec[a.set]  = p.ncpen.fun(y.vec,ax.mat,b.vec[a.set],w.vec[a.set],lam,gam,tau,iter.max,b.eps,k.eps,p.eff,r.eff)$b.vec
                    // cpp
                    //b.vec[a.set] = p.ncpen.fun(y.vec,ax.mat,b.vec[a.set],w.vec[a.set],lam,gam,tau,alp,iter.max,b.eps,k.eps,p.eff,fam,pen)$b.vec
                    arma::mat ax_mat;
                    if(pen == "ridge") {
                         p_ncpen_ret ret_buff;
                         p_ncpen_fun(y_vec, x_mat, b_vec, w_vec, lam, gam, tau, alp, iter_max, qiter_max, qiiter_max,
                                     b_eps, k_eps, p_eff, cut, c_eps,
                                     obj_fun, obj_grad_fun, obj_hess_fun, pen_fun, pen_grad_fun, ret_buff);
                         b_vec = ret_buff.b_vec;
                    } else {
                         if(fam == "cox") {
                              arma::uvec ca_set = a_set;
                              ca_set.resize(ca_set.n_rows+1);
                              // index starts from 0.
                              // At the beginning, p = p-1.
                              // So, p is the index for the last columin of x_mat which contains censored indicators.
                              ca_set(ca_set.n_rows-1) = p;
                              ax_mat = x_mat.cols(ca_set);
                         } else {
                              ax_mat = x_mat.cols(a_set);
                         }

                         arma::vec ab_vec = b_vec(a_set);
                         arma::vec aw_vec = w_vec(a_set);
                         // b_vec(a_set) = Rcpp::as<arma::vec>(
                         //      p_ncpen_fun(y_vec, ax_mat, ab_vec, aw_vec, lam, gam, tau, iter_max, b_eps, k_eps, p_eff, r_eff,
                         //                  obj_fun, obj_grad_fun, obj_hess_fun, pen_fun, pen_grad_fun)["b.vec"]);
                         //b_vec(a_set).print();

                         p_ncpen_ret ret_buff;
                         p_ncpen_fun(y_vec, ax_mat, ab_vec, aw_vec, lam, gam, tau, alp, iter_max, qiter_max, qiiter_max,
                                     b_eps, k_eps, p_eff, cut, c_eps,
                                     obj_fun, obj_grad_fun, obj_hess_fun, pen_fun, pen_grad_fun, ret_buff);
                         b_vec(a_set) = ret_buff.b_vec;
                    }
                    //#############################################################################################################

                    //g.vec = obj.grad.fun(y.vec,x.mat,b.vec) + alp*w.vec*pen.grad.fun(b.vec,lam,gam,tau) - alp*lam*w.vec*sign(b.vec) + 2*(1-alp)*lam*w.vec*b.vec
                    g_vec = obj_grad_fun(y_vec, x_mat, b_vec) + alp*w_vec%pen_grad_fun(b_vec, lam, gam, tau) - alp*lam*w_vec%arma::sign(b_vec) + 2*(1-alp)*lam*w_vec%b_vec;


                    //p.vec = alp*lam*w.vec*sign(b.vec);
                    p_vec = alp*lam*w_vec%arma::sign(b_vec);
                    //a.set = c(1:p)[b.vec!=0]
                    a_set = arma::find(b_vec!=0);
                    //n.set = c(1:p)[b.vec==0]
                    n_set = arma::find(b_vec==0);
                    //kkt0 = sum(abs(g.vec[a.set]+p.vec[a.set])<k.eps)==length(a.set)     # Mod sum(abs( )<k.eps)==  a.set; last sum => length
                    kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps) == a_set.n_rows;
                    //kkt1 = sum(abs(g.vec[n.set])-alp*lam*w.vec[n.set]<k.eps)==length(n.set)
                    kkt1 = arma::sum(arma::abs(g_vec(n_set))-alp*lam*w_vec(n_set)<k_eps)== n_set.n_rows;

                    if(NCPEN_DEVELOP_MODE == true) {
                         //bool con2 = kkt0&&kkt1;
                         Rcpp::Rcout << "______________________________________________________ncpen_fun : a_set=" << a_set.n_rows <<
                              ", lambda=" << lam << ", iter=" << iter << ", kkt0=" << kkt0 << ", kkt1=" << kkt1 << std::endl;
                    }
                    //if(kkt0&kkt1) break #iter break
                    if(kkt0&&kkt1) { break; } // iter break
               } //iter
          }

          if(NCPEN_DEVELOP_MODE == true) {
               //bool con2 = kkt0&&kkt1;
               Rcpp::Rcout << "______________________________________________________ncpen_fun : a_set=" << a_set.n_rows <<
                    ", lambda=" << lam << ", pos=" << pos << ", kkt0=" << kkt0 << ", kkt1=" << kkt1 << std::endl;
          }

          //b.mat[,pos] = b.vec
          b_mat.col(pos) = b_vec;
          //g.mat[,pos] = g.vec+p.vec
          g_mat.col(pos) = g_vec + p_vec;
          //c.mat[,pos] = c(kkt0,kkt1)
          c_mat(0, pos) = kkt0; c_mat(1, pos) = kkt1;
          //f.vec[pos] = obj.fun(y.vec,x.mat,b.vec) + alp*sum(w.vec*pen.fun(b.vec,lam,gam,tau)) + (1-alp)*lam*sum(w.vec*b.vec^2)
          f_vec(pos) = obj_fun(y_vec, x_mat, b_vec) + alp*arma::sum(w_vec%pen_fun(b_vec, lam, gam, tau)) + (1-alp)*lam*arma::sum(w_vec%arma::pow(b_vec,2));

          //d.vec[pos] = length(a.set)    ### (Mod) sum => length
          a_set = arma::find(b_vec!=0);
          d_vec(pos) = a_set.n_rows;

          // Mod
          //if(length(a.set)>d.max){
          if(a_set.n_rows > d_max) {
               //cat("(warning) ncpen stops because the number of non-zero parameters exceeds d.max", d.max,"\n")
               //Rcout << "(Warning) ncpen_fun stops because the number of non-zero parameters exceeds d_max " << d_max << "." << std::endl;
               //cat("(warning) increase d.max", "\n")                                        ### (Mod)
               //Rcout << "(Warning) increase d_max." << std::endl;
               warning_code_vec(0) = 1;
               //break #pos break
               pos++; // for head_cols functions
               break; // pos break
          }

          // Add
          // if(fam=="log"){
          if(obj_fun == log_obj_fun) {
               //if(max(abs(b.vec))>div){
               if(arma::max(arma::abs(b_vec))>div){
                    // cat("(warning) ncpen stops because a paramter becomes larger than div=",div,"\n")
                    // cat("(warning) try larger ridge effect by reducing current alp=",alp, "\n")
                    //Rcout << "(Warning) ncpen_fun stops because a paramter becomes larger than div=" << div << std::endl;
                    //Rcout << "(Warning) try larger ridge effect by reducing current alp=" << alp << std::endl;
                    warning_code_vec(1) = 1;
                    //break #pos break
                    pos++; // for head_cols functions
                    break; //pos break
               }
          }
          // Add

     } // pos
     // return(list(b.mat=b.mat[,1:pos],g.mat=g.mat[,1:pos],f.vec=f.vec[1:pos],c.mat=c.mat[,1:pos],
     //             lam.vec=lam.vec[1:pos],d.vec=d.vec[1:pos],w.vec=w.vec))                         ### (Mod) add w.vec in return

     arma::mat ret_b_mat = b_mat.head_cols(pos);
     arma::mat ret_g_mat = g_mat.head_cols(pos);
     arma::vec ret_f_vec = f_vec.head(pos);
     arma::mat ret_c_mat = c_mat.head_cols(pos);
     arma::vec ret_lam_vec = lam_vec.head(pos);
     arma::vec ret_d_vec = d_vec.head(pos);

     // return Rcpp::List::create(Rcpp::Named("beta") = ret_b_mat,
     //                           Rcpp::Named("gradient") = ret_g_mat,
     //                           Rcpp::Named("f.vec") = ret_f_vec,
     //                           Rcpp::Named("conv") = ret_c_mat,
     //                           Rcpp::Named("lambda") = ret_lam_vec,
     //                           Rcpp::Named("df") = ret_d_vec,
     //                           Rcpp::Named("w.lambda") = w_vec,
     //                           Rcpp::Named("warnings") = warning_code_vec);

     ret_buff.b_mat = ret_b_mat;
     ret_buff.g_mat = ret_g_mat;
     ret_buff.f_vec = ret_f_vec;
     ret_buff.c_mat = ret_c_mat;
     ret_buff.lam_vec = ret_lam_vec;
     ret_buff.d_vec = ret_d_vec;
     ret_buff.w_vec = w_vec;
     ret_buff.warnings = warning_code_vec;

     // std::cout << "CPP_NCPEN ver 0.1.9.0" << std::endl;

     return 0;

}

/* End of main functions --------------------------------------------------------*/
