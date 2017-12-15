/*--------------------------------------------------------
 * Includes
 --------------------------------------------------------*/
#include "ncpen.h"
#include <RcppArmadillo.h>

using namespace arma;

/* End of includes --------------------------------------------------------*/

/*--------------------------------------------------------
 * Utility functions
 --------------------------------------------------------*/
template <typename T> int sign(T val) {
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

arma::vec rm_row(arma::uword rm_idx, arma::vec& v) {
     //     Rcout << "rm_idx: " << rm_idx << std::endl;
     arma::uword len = v.n_rows;
     arma::vec ret(len-1);
     ret.head(rm_idx) = v.head(rm_idx);
     ret.tail(len - rm_idx - 1) = v.tail(len - rm_idx - 1);

     return ret;
}


arma::vec rm_row(arma::vec rm_vec, arma::uword start, arma::uword end) {
     rm_vec = sort(rm_vec);
     arma::uword len = (end-start + 1) - rm_vec.n_rows;
     arma::vec ret(len);

     arma::uword ii = 0;
     arma::uword rm_i = 0;
     for(arma::uword i = 0; i<len; i++) {
          if(rm_vec(rm_i) != i) {
               ret(ii) = i;
               ii++;
          } else {
               rm_i++;
          }
     }

     return ret;
}


int which_max(arma::vec& v) { //do not use uword to return -1
     if(v.n_rows <=0) return -1;

     arma::uword max_idx = 0;
     double max_val = v[0];

     for(arma::uword i=1; i<v.n_rows; i++) {
          if(v(i) > max_val) {
               max_val = v(i);
               max_idx = i;
          }
     }

     return max_idx;
}

arma::uvec rm_row(arma::uvec rm_vec, arma::uword start, arma::uword end) {
     rm_vec = sort(rm_vec);
     arma::uword len = (end-start + 1) - rm_vec.n_rows;
     arma::uvec ret(len);

     arma::uword ii = 0;
     arma::uword rm_i = 0;
     for(arma::uword i = 0; ; i++) {
          //Rcout << "i: " <<i<<"  ii: " <<ii<<"  rm_i: " <<rm_i <<endl;
          if(rm_i == rm_vec.n_rows) {
               for(;ii<len;) {
                    ret(ii) = i;
                    ii++;
                    i++;
               }
               break;
          }

          if(rm_vec(rm_i) != i) {
               ret(ii) = i;
               ii++;
          } else {
               rm_i++;
          }
     }

     return ret;
}

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
//
// if(pen=="sridge"){
//      pen.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec)
//           tem0 = (lam-2*gam*tau*(lam-gam))/(tau^2*(lam-gam)^2)
//           tem1 = -tem0*tau^3*(lam-gam)^3/3+lam*tau*(lam-gam)-gam*tau^2*(lam-gam)^2
//           tem0.vec = (-tem0*ab.vec^3/3+lam*ab.vec) * (ab.vec<tau*(lam-gam))
//           tem1.vec = (gam*ab.vec^2+tem1) * (ab.vec>=tau*(lam-gam))
//           return(tem0.vec+tem1.vec)
//      }
//      pen.grad.fun = function(b.vec,lam,gam,tau){
//           ab.vec = abs(b.vec); sb.vec = sign(b.vec)
//                tem0 = (lam-2*gam*tau*(lam-gam))/(tau^2*(lam-gam)^2)
//                tem0.vec = (-tem0*ab.vec^2+lam) * (ab.vec<tau*(lam-gam))
//                tem1.vec = (2*gam*ab.vec) * (ab.vec>=tau*(lam-gam))
//                return((tem0.vec+tem1.vec)*sb.vec)
//      }
// }
arma::vec sridge_pen_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     double tem0 = (lam-2*gam*tau*(lam-gam))/(std::pow(tau,2)*pow(lam-gam,2));
     double tem1 = -tem0*std::pow(tau,3)*pow(lam-gam,3)/3+lam*tau*(lam-gam)-gam*std::pow(tau,2)*std::pow(lam-gam,2);
     arma::vec tem0_vec = (-tem0*arma::pow(ab_vec,3)/3+lam*ab_vec) % (ab_vec<tau*(lam-gam));
     arma::vec tem1_vec = (gam*arma::pow(ab_vec,2)+tem1) % (ab_vec>=tau*(lam-gam));

     return tem0_vec + tem1_vec;
}
arma::vec sridge_pen_grad_fun(arma::vec& b_vec, double lam, double gam, double tau) {
     arma::vec ab_vec = arma::abs(b_vec);
     //     arma::vec sb_vec = arma::sign(b_vec);
     double tem0 = (lam-2*gam*tau*(lam-gam))/(std::pow(tau,2)*std::pow(lam-gam,2));
     arma::vec tem0_vec = (-tem0*arma::pow(ab_vec,2)+lam) % (ab_vec<tau*(lam-gam));
     arma::vec tem1_vec = (2*gam*ab_vec) % (ab_vec>=tau*(lam-gam));

     //     return (tem0_vec+tem1_vec)%sb_vec;
     return (tem0_vec+tem1_vec)%arma::sign(b_vec);
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
     } else {
          throw std::invalid_argument("Invalid penalty gradient funtion option. Only available \"scad\", \"mcp\", \"tlp\", \"classo\", \"sridge\", \"mbridge\", \"mlog\" or \"lasso\".");
          return NULL;
     }
}
// End of penalty functions ------------------------------------------


// 2. Loss functions ------------------------------------------

// 1. if(fam=="lin"){
double lin_obj_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = drop(x.mat%*%b.vec)
     arma::vec xb_vec = x_mat*b_vec;
//   return(  sum((xb.vec-y.vec)^2)/length(y.vec)/2+ r.eff*sum(b.vec^2) )    #(Modified) -sum=>sum; y-xb=>xb-y; r.eff
     return arma::sum(arma::pow(xb_vec-y_vec, 2))/y_vec.n_rows/2 + r_eff*arma::sum(b_vec%b_vec); //(Modified) -sum=>sum; y-xb=>xb-y; r.eff
}

arma::vec lin_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = drop(x.mat%*%b.vec)
     arma::vec xb_vec = x_mat*b_vec;
//   return( drop(t(x.mat)%*%(xb.vec-y.vec)/length(y.vec))+ 2*r.eff*b.vec )  #(Modified) -t=>t; y-xb=>xb-y; r.eff
     return x_mat.t()*(xb_vec-y_vec)/y_vec.n_rows+ 2*r_eff*b_vec;  //(Modified) -t=>t; y-xb=>xb-y; r.eff
}

//??
arma::mat lin_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   d.mat = diag(rep(1,length(b.vec)))                         #(Added) d.mat= diagonal matrix
     arma::mat d_mat = arma::eye<mat>(b_vec.n_rows, b_vec.n_rows);
//   return( t(x.mat)%*%x.mat/length(y.vec) + 2*r.eff*d.mat )   #(Mod) r.eff; d.mat
     return x_mat.t()*x_mat/y_vec.n_rows + 2*r_eff*d_mat; //(Mod) r.eff; d.mat
}


// 2. if(fam=="poi"){

double poi_obj_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = pmin(drop(x.mat%*%b.vec),700);                                      #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>POI_BIG_XB)).fill(POI_BIG_XB);
//   return( sum(exp(xb.vec)-y.vec*xb.vec)/length(y.vec) + 2*r.eff*sum(b.vec^2) ) #(Mod) r.eff
     return arma::sum(arma::exp(xb_vec)-y_vec%xb_vec)/y_vec.n_rows + 2*r_eff*arma::sum(b_vec%b_vec); //(Mod) r.eff
}

arma::vec poi_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = pmin(drop(x.mat%*%b.vec),700);                                      #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>POI_BIG_XB)).fill(POI_BIG_XB);
     //   return( drop(t(x.mat)%*%(exp(xb.vec)-y.vec))/length(y.vec) + 2*r.eff*b.vec ) #(Mod) -drop=>drop; y-xb=>exp(xb)-y; r.eff
     return x_mat.t()*(arma::exp(xb_vec)-y_vec)/y_vec.n_rows + 2*r_eff*b_vec; //(Mod) -drop=>drop; y-xb=>exp(xb)-y; r.eff
}

arma::mat poi_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = pmin(drop(x.mat%*%b.vec),700);                                   # Mod 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>POI_BIG_XB)).fill(POI_BIG_XB);
     //   exb.vec = exp(xb.vec)                                                     # Add exb.vec
     arma::vec exb_vec = arma::exp(xb_vec); // Add exb.vec
     //   d.mat = diag(rep(1,length(b.vec)))                                        # Add d.mat=´ë°¢¿ø¼Ò ¸ðµÎ 1 ³ª¸ÓÁö 0ÀÎ ¸ÅÆ®¸¯½º
     arma::mat d_mat = arma::eye<mat>(b_vec.n_rows, b_vec.n_rows);
//   return( t(x.mat)%*%diag(exb.vec)%*%x.mat/length(y.vec) + 2*r.eff*d.mat )  # Mod exb.vec; r.eff; d.mat
     return x_mat.t()*arma::diagmat(exb_vec)*x_mat/y_vec.n_rows + 2*r_eff*d_mat; // Mod exb.vec; r.eff; d.mat
}

// 3. if(fam=="log"){
double log_obj_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = pmin(drop(x.mat%*%b.vec),700)                                              # Mod 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>LOGI_BIG_XB)).fill(LOGI_BIG_XB);
     //   return( sum(log(1+exp(xb.vec))-y.vec*xb.vec)/length(y.vec) + r.eff*sum(b.vec^2) )   # Mod -sum=>sum; y-log=>log-y; r.eff
     return arma::sum(arma::log(1+arma::exp(xb_vec))-y_vec%xb_vec)/y_vec.n_rows + r_eff*arma::sum(b_vec%b_vec); // Mod -sum=>sum; y-log=>log-y; r.eff
}

arma::vec log_obj_grad_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
//   xb.vec = pmin(drop(x.mat%*%b.vec),700)                                  #(Mod) 300 => 700
     arma::vec xb_vec = x_mat*b_vec;
     xb_vec(arma::find(xb_vec>LOGI_BIG_XB)).fill(LOGI_BIG_XB);
//   exb.vec = exp(xb.vec)
     arma::vec exb_vec = arma::exp(xb_vec);
//   p.vec = exb.vec/(1+exb.vec)
     arma::vec p_vec = exb_vec/(1+exb_vec);
//   return( drop(t(x.mat)%*%(p.vec-y.vec)/length(y.vec)) + 2*r.eff*b.vec)   #(Mod) r.eff
     return x_mat.t()*(p_vec-y_vec)/y_vec.n_rows + 2*r_eff*b_vec;
}

arma::mat log_obj_hess_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec, double r_eff) {
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
//   d.mat = diag(rep(1,length(b.vec)))                                       #(Add) d.mat=´ë°¢¿ø¼Ò ¸ðµÎ 1 ³ª¸ÓÁö 0ÀÎ ¸ÅÆ®¸¯½º
     arma::mat d_mat = arma::eye<mat>(b_vec.n_rows, b_vec.n_rows);
//   return( t(x.mat)%*%diag(d.vec)%*%x.mat/length(y.vec) + 2*r.eff*d.mat )   #(Mod) exb.vec; r.eff; d.mat
     return x_mat.t()*arma::diagmat(d_vec)*x_mat/y_vec.n_rows + 2*r_eff*d_mat; //(Mod) exb.vec; r.eff; d.mat
}

obj_fun_ptr get_obj_fun_ptr(std::string name) {
     if(name.compare("gaussian") == 0) {
          return lin_obj_fun;
     } else if(name.compare("poisson") == 0) {
          return poi_obj_fun;
     } else if(name.compare("binomial") == 0) {
          return log_obj_fun;
     } else {
          throw std::invalid_argument("Invalid oject function family. Only supports \"gaussian\" (linear), \"binomial\" (logistic), \"poisson\".");
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
     } else {
          throw std::invalid_argument("Invalid oject function family @ gradient. Only supports \"gaussian\" (linear), \"binomial\" (logistic), \"poisson\".");
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
     } else {
          throw std::invalid_argument("Invalid oject function family @ hessian. Only supports \"gaussian\" (linear), \"binomial\" (logistic), \"poisson\".");
          return NULL;
     }
}
// End of loss functions ------------------------------------------
// ###### Newton-Raphosn for nun-penalized model : 2017.0803
// nr.fun = function(y.vec,x.mat,r.eff,iter.max,b.eps){
arma::vec nr_fun(std::string fam, arma::vec& y_vec, arma::mat& x_mat, double r_eff, double iter_max, double b_eps) {
     //std::string fam = Rcpp::as<std::string>(family);

     //   b.vec = rep(0,dim(x.mat)[2])
     arma::vec b_vec = arma::zeros<arma::vec>(x_mat.n_cols);
     //   for(iter in 1:iter.max){
     arma::uword iter = 0;
     for(;iter <iter_max; iter++) {
          // #print(b.vec)
          //        hess.mat = obj.hess.fun(y.vec,x.mat,b.vec,r.eff)
          arma::mat hess_mat = get_obj_hess_fun_ptr(fam)(y_vec, x_mat, b_vec, r_eff);
          //        grad.vec = obj.grad.fun(y.vec,x.mat,b.vec,r.eff)
          arma::vec grad_vec = get_obj_grad_fun_ptr(fam)(y_vec, x_mat, b_vec, r_eff);
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
     return sign(est)*(std::abs(est)-del)*(std::abs(est)>del);
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
     arma::vec q_mat_col = q_mat.col(pos);
     return -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;
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
                                  double lam, double iter_max, double b_eps, double k_eps, arma::uword p_eff, arma::uword q_rank, p_ncpen_ret& ret_buff) {

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
     // Vriables for loop resue----------------------

     arma::uword iter = 0;
//   for(iter in 1:iter.max){#iter
     for(;iter < iter_max; iter++){ // iter
// ##### active set
          //p0.eff = p.eff                                  ### (Add)
          arma::uword pp_eff = p_eff; //### p0.eff => pp.eff
          //ob.vec = b.vec                                  ### (Mod) oob.vec => ob.vec
          ob_vec = b_vec; // (Mod) oob.vec => ob.vec
          //a1.set = c(1:p)[(b.vec!=0)&(w.vec!=0)]
          arma::uvec a1_set = arma::find( ((b_vec!= 0) % (w_vec!=0))> 0 );
          //a2.set = c(1:p)[(b.vec!=0)&(w.vec==0)]
          arma::uvec a2_set = arma::find( ((b_vec!= 0) % (w_vec==0))> 0 );
          //for(iiter in 1:iter.max){#iiter
          //for(arma::uword iiter = 0; iiter < iter_max; iiter++) { // iiter
          arma::uword iiter = 0;
          for(; iiter < iter_max; iiter++) { // iiter
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
               //if(sum(abs(b.vec-ob.vec))<b.eps) break           ### (Mod) oob.vec => ob.vec
               if(arma::sum(arma::abs(b_vec-ob_vec))<b_eps) {break;}
               //#####################################
               //##### projection: reduced model #####
               //#####################################
               //b.vec[abs(b.vec)<b.eps] = 0                                ### add this
               b_vec(arma::find(arma::abs(b_vec) < b_eps)).fill(0);
               //ob.vec = b.vec                                   ### (Mod) oob.vec => ob.vec
               ob_vec = b_vec;

               //############################################# Add
               //##### iiter projection
               //if(iiter>p0.eff){
               if(iiter > pp_eff){ //### p0.eff => pp.eff
                    //a.set = c(1:p)[b.vec!=0]
                    arma::uvec a_set = arma::find(b_vec != 0);
                    //if(length(a.set)<q.rank){
                    if(a_set.n_rows < q_rank){ //### please check the inequality here! if length(a.set)>q.rank=n then the inverse may not exist!!!
                         //pb.vec = b.vec*0
                         arma::vec pb_vec = arma::zeros<arma::vec>(b_vec.n_rows);

                         //[s 2017011]
                         // //pb.vec[a.set] = -solve(q.mat[a.set,a.set])%*%(l.vec[a.set]+lam*w.vec[a.set]*sign(b.vec[a.set]))/2
                         // pb_vec(a_set) = -arma::inv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                         //
                         // //pbf = t(pb.vec)%*%(q.mat)%*%pb.vec + sum(l.vec*pb.vec)+lam*sum(abs(w.vec*pb.vec))
                         // double pbf = arma::as_scalar(pb_vec.t()*q_mat*pb_vec) + arma::sum(l_vec%pb_vec)+lam*arma::sum(arma::abs(w_vec%pb_vec));
                         //
                         // //bf  = t( b.vec)%*%(q.mat)%*% b.vec + sum(l.vec* b.vec)+lam*sum(abs(w.vec* b.vec))
                         // double bf = arma::as_scalar(b_vec.t()*q_mat*b_vec) + arma::sum(l_vec%b_vec)+lam*arma::sum(arma::abs(w_vec%b_vec));

                         //pb.vec[a.set] = -ginv(q.mat[a.set,a.set])%*%(l.vec[a.set]+lam*w.vec[a.set]*sign(b.vec[a.set]))/2                  #### solve -> ginv
                         //pb_vec(a_set) = -arma::pinv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                         pb_vec(a_set) = -arma::inv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;//pbf = t(pb.vec[a.set])%*%(q.mat[a.set,a.set])%*%pb.vec[a.set] + sum(l.vec[a.set]*pb.vec[a.set])+lam*sum(abs(w.vec[a.set]*pb.vec[a.set]))
                         double pbf = arma::as_scalar(pb_vec(a_set).t()*q_mat(a_set, a_set)*pb_vec(a_set)) + arma::sum(l_vec(a_set)%pb_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%pb_vec(a_set)));
                         //bf  = t( b.vec[a.set])%*%(q.mat[a.set,a.set])%*% b.vec[a.set] + sum(l.vec[a.set]* b.vec[a.set])+lam*sum(abs(w.vec[a.set]* b.vec[a.set]))
                         double bf = arma::as_scalar(b_vec(a_set).t()*q_mat(a_set, a_set)*b_vec(a_set)) + arma::sum(l_vec(a_set)%b_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%b_vec(a_set)));
                         //[e 2017011]

                         //if(pbf<bf){ ob.vec = b.vec = pb.vec; cat("iiter projenction works","\n")
                         if(pbf < bf){ // How about if((pbf-bf)<b_eps) {
                              //ob.vec = b.vec = pb.vec;
                              ob_vec = pb_vec;
                              b_vec = pb_vec;
                              //pp.eff = p.eff                                                      ### add
                              pp_eff = p_eff;                                                      //### add
                              //cat("iiter projenction works","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                              //Rcout << "iiter projection worked. lam=" << lam << ", df=" << sum(b_vec!=0) << std::endl;
                         //} else { p0.eff = 2*p0.eff; cat("iiter projenction fails","\n");  }
                         } else {
                              //p0.eff = 2*p0.eff;
                              pp_eff = 2*pp_eff;
                              //cat("iiter projenction fails","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                              //Rcout << "iiter projection failed. lam" << lam << ", df=" << sum(b_vec!=0) << std::endl;
                         }
                    }
               }
               //############################################# Add

          }//iiter
          //Rcout << "CPP native_qlasso_fun | iiter: " << iiter << std::endl;

          //b.vec[abs(b.vec)<1e-7] = 0      ### (Add)
          //[s/e 20170111] b_vec(arma::find(arma::abs(b_vec) < SMALL_NUMBER)).fill(0);

          //###################
          //##### null set
          //###################
          //n1.set = c(1:p)[(b.vec==0)&(w.vec!=0)]
          arma::uvec n1_set = arma::find( ((b_vec== 0) % (w_vec!=0))> 0 );
          //n2.set = c(1:p)[(b.vec==0)&(w.vec==0)]
          arma::uvec n2_set = arma::find( ((b_vec== 0) % (w_vec==0))> 0 );
          //for(pos in n1.set){
          for(arma::uword idx = 0; idx < n1_set.n_rows; idx++) {
               arma::uword pos = n1_set(idx);
               //est = -(2*sum(q.mat[-pos,pos]*b.vec[-pos])+l.vec[pos])/q.mat[pos,pos]/2
               //double est = -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;
               double est = get_qlasso_fun_est(pos, q_mat, b_vec, l_vec);
               //del = lam*w.vec[pos]/q.mat[pos,pos]/2
               //double del = lam*w_vec(pos)/q_mat(pos, pos)/2;
               double del = get_qlasso_fun_del(pos, q_mat, w_vec, lam);
               //b.vec[pos] = soft.fun(est,del)
               b_vec(pos) = soft_fun(est, del);
          }

          //for(pos in n2.set){
          for(arma::uword idx = 0; idx < n2_set.n_rows; idx++) {
               arma::uword pos = n2_set(idx);
               //est = -(2*sum(q.mat[-pos,pos]*b.vec[-pos])+l.vec[pos])/q.mat[pos,pos]/2
               //double est = -(2*arma::sum(rm_row(pos, q_mat_col)%rm_row(pos, b_vec))+l_vec(pos))/q_mat(pos, pos)/2;
               double est = get_qlasso_fun_est(pos, q_mat, b_vec, l_vec);
               //b.vec[pos] = est
               b_vec(pos) = est;
          }

          //b.vec[abs(b.vec)<1e-7] = 0    ### (Add)
          b_vec(arma::find(arma::abs(b_vec) < b_eps)).fill(0); //### 1e-7 => b.eps


          // ##################################
          // ##### projection: full model #####
          // ##################################
          //if(iter>p.eff){
          if(iter > p_eff){ //### please check the inequality here! if length(a.set)>q.rank=n then the inverse may not exist!!!
               //a.set = c(1:p)[b.vec!=0]
               arma::uvec a_set = arma::find(b_vec != 0);
               //if(length(a.set)<q.rank){
               if(a_set.n_rows<q_rank){
                    //pb.vec = b.vec*0
                    arma::vec pb_vec = b_vec*0;
                    //[s 20170111]
                    // //pb.vec[a.set] = -solve(q.mat[a.set,a.set])%*%(l.vec[a.set]+lam*w.vec[a.set]*sign(b.vec[a.set]))/2        ### (Mod) Check parentheses
                    // pb_vec(a_set) = -arma::inv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                    // //pbf = t(pb.vec)%*%(q.mat)%*%pb.vec + sum(l.vec*pb.vec)+lam*sum(abs(w.vec*pb.vec))
                    // double pbf = arma::as_scalar(pb_vec.t()*q_mat*pb_vec) + arma::sum(l_vec%pb_vec)+lam*arma::sum(arma::abs(w_vec%pb_vec));
                    // //bf  = t( b.vec)%*%(q.mat)%*% b.vec + sum(l.vec* b.vec)+lam*sum(abs(w.vec* b.vec))
                    // double bf = arma::as_scalar(b_vec.t()*q_mat*b_vec) + arma::sum(l_vec%b_vec)+lam*arma::sum(arma::abs(w_vec%b_vec));

                    //pb.vec[a.set] = -ginv(q.mat[a.set,a.set])%*%(l.vec[a.set]+lam*w.vec[a.set]*sign(b.vec[a.set]))/2                  ### solve => ginv
                    //pb_vec(a_set) = -arma::pinv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                    pb_vec(a_set) = -arma::inv(q_mat(a_set, a_set))*(l_vec(a_set)+lam*w_vec(a_set)%arma::sign(b_vec(a_set)))/2;
                    ///pbf = t(pb.vec[a.set])%*%(q.mat[a.set,a.set])%*%pb.vec[a.set] + sum(l.vec[a.set]*pb.vec[a.set])+lam*sum(abs(w.vec[a.set]*pb.vec[a.set]))
                    double pbf = arma::as_scalar(pb_vec(a_set).t()*q_mat(a_set, a_set)*pb_vec(a_set)) + arma::sum(l_vec(a_set)%pb_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%pb_vec(a_set)));
                    //bf  = t( b.vec[a.set])%*%(q.mat[a.set,a.set])%*% b.vec[a.set] + sum(l.vec[a.set]* b.vec[a.set])+lam*sum(abs(w.vec[a.set]* b.vec[a.set]))
                    double bf = arma::as_scalar(b_vec(a_set).t()*q_mat(a_set, a_set)*b_vec(a_set)) + arma::sum(l_vec(a_set)%b_vec(a_set))+lam*arma::sum(arma::abs(w_vec(a_set)%b_vec(a_set)));
                    //[e 20170111]

                    if(pbf<bf){
                         //ob.vec = b.vec = pb.vec;
                         ob_vec = b_vec = pb_vec;
                         //cat("iter projenction works","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                         //Rcout << "iter projection worked. lam=" << lam << ", df=" << sum(b_vec!=0) << std::endl;
                    } else {
                         //p.eff = 2*p.eff;
                         p_eff = 2*p_eff;
                         //cat("iter projenction fails","lam=",lam,"df=",sum(b.vec!=0),"\n")  ### add "lam="
                         //Rcout << "iter projection failed. lam" << lam << ", df=" << sum(b_vec!=0) << std::endl;
                    }
               }
          }
          //##### iter projection

          //b.vec[abs(b.vec)<1e-7] = 0 ### (Add)
          //[s/e 20170111] b_vec(arma::find(arma::abs(b_vec) < SMALL_NUMBER)).fill(0);

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
          //a.set = c(1:p)[b.vec!=0]                                               #### (Add)
          arma::uvec a_set = arma::find(b_vec!=0); //                                               #### (Add
          //n.set = c(1:p)[b.vec==0]                                               #### (Add)
          arma::uvec n_set = arma::find(b_vec==0); //                                               #### (Add
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
//p.ncpen.fun = function(y.vec,x.mat,b.vec,w.vec,lam,gam,tau,iter.max,b.eps,k.eps,p.eff,r.eff){ ## (Mod) r.eff
// Rcpp::List p_ncpen_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec0, arma::vec& w_vec,
//                        double lam, double gam, double tau, double iter_max, double b_eps, double k_eps, arma::uword p_eff, double r_eff,
//                        obj_fun_ptr obj_fun, obj_grad_fun_ptr obj_grad_fun, obj_hess_fun_ptr obj_hess_fun,
//                        pen_fun_ptr pen_fun, pen_grad_fun_ptr pen_grad_fun) {
int p_ncpen_fun(arma::vec& y_vec, arma::mat& x_mat, arma::vec& b_vec0, arma::vec& w_vec,
                       double lam, double gam, double tau, double iter_max, double b_eps, double k_eps, arma::uword p_eff, double r_eff,
                       obj_fun_ptr obj_fun, obj_grad_fun_ptr obj_grad_fun, obj_hess_fun_ptr obj_hess_fun,
                       pen_fun_ptr pen_fun, pen_grad_fun_ptr pen_grad_fun, p_ncpen_ret& ret_buff) {

     //q.rank = length(y.vec)        ### (Mod)
     int q_rank = y_vec.n_rows;
     //p = length(b.vec)             ### (Add)
     arma::vec b_vec = b_vec0;
     //int p = b_vec.n_rows;
     //ob.vec = b.vec
     arma::vec ob_vec = b_vec;

     arma::vec g_vec(b_vec.n_rows);
     arma::vec p_vec(b_vec.n_rows);
     bool kkt0 = false;
     bool kkt1 = false;
     //f.vec = rep(0,iter.max)
     arma::vec f_vec = zeros<vec>(iter_max);

     //for(iter in 1:iter.max){
     arma::uword iter = 0;
     for(; iter<iter_max; iter++) {
          //h.mat = obj.hess.fun(y.vec,x.mat,b.vec,r.eff)      ## (Mod) r.eff
          arma::mat h_mat = obj_hess_fun(y_vec, x_mat, b_vec, r_eff);
          //g.vec = obj.grad.fun(y.vec,x.mat,b.vec,r.eff)      ## (Mod) r.eff
          g_vec = obj_grad_fun(y_vec, x_mat, b_vec, r_eff);
          //l.vec = g.vec - drop(h.mat%*%b.vec)+w.vec*pen.grad.fun(b.vec,lam,gam,tau)-lam*w.vec*sign(b.vec)   ## (Mod)
          arma::vec l_vec = g_vec - h_mat*b_vec + w_vec%pen_grad_fun(b_vec, lam, gam, tau)-lam*w_vec%arma::sign(b_vec);
          //b.vec = qlasso.fun(h.mat/2,l.vec,b.vec,w.vec,lam,iter.max,b.eps,k.eps,p.eff,q.rank)$b.vec
          arma::mat h_mat_2 = h_mat/2;
          //b_vec = Rcpp::as<arma::vec>(qlasso_fun(h_mat_2, l_vec, b_vec, w_vec, lam, iter_max, b_eps, k_eps, p_eff, q_rank)["b.vec"]);
          p_ncpen_ret ret_buff;
          qlasso_fun(h_mat_2, l_vec, b_vec, w_vec, lam, iter_max, b_eps, k_eps, p_eff, q_rank, ret_buff);
          b_vec = ret_buff.b_vec;

          //t.vec = w.vec*(pen.grad.fun(ob.vec,lam,gam,tau)-lam*sign(ob.vec))                          ### (Add)  use ob.vec
          arma::vec t_vec = w_vec%(pen_grad_fun(ob_vec, lam, gam, tau)-lam*arma::sign(ob_vec)); // (Add)  use ob.vec
          //oob = obj.fun(y.vec,x.mat,ob.vec,r.eff) + sum(t.vec*ob.vec) + lam*sum(w.vec*abs(ob.vec))   ### (Add)  use ob.vec
          double oob = obj_fun(y_vec, x_mat, ob_vec, r_eff) + arma::sum(t_vec%ob_vec) + lam*arma::sum(w_vec%arma::abs(ob_vec)); // (Add)  use ob.vec
          //ob  = obj.fun(y.vec,x.mat, b.vec,r.eff) + sum(t.vec* b.vec) + lam*sum(w.vec*abs( b.vec))   ### (Add)  use b.vec
          double ob  = obj_fun(y_vec, x_mat,  b_vec, r_eff) + arma::sum(t_vec% b_vec) + lam*arma::sum(w_vec%arma::abs( b_vec)); // (Add)  use b.vec

          //#################################################### if - for Add
          //#################################################### if - for Add
          //if(ob>oob){#mlqa
          //if(ob>oob) { // mlqa
          //if(ob>(oob+b.eps)){#mlqa
          if(ob>oob+b_eps) {
               //cat("start golden section algorithm","\n")
               //Rcout << "starting golden section algorithm..." << std::endl;
               //for(iiter in 1:iter.max){
               //for(arma::uword iiter = 0; iiter < iter_max; iiter++) { // iiter
               arma::uword iiter = 0;
               for(; iiter < iter_max; iiter++) { // iiter
                    //gold = (sqrt(5)-1)/2
                    double gold = (std::sqrt(5)-1)/2;
                    //b1.vec = (1-gold)* b.vec + gold*ob.vec
                    arma::vec b1_vec = (1-gold)* b_vec + gold*ob_vec;
                    //b2.vec = (1-gold)*ob.vec + gold* b.vec
                    arma::vec b2_vec = (1-gold)*ob_vec + gold* b_vec;
                    //ob1 = obj.fun(y.vec,x.mat,b1.vec,r.eff) + sum(t.vec*b1.vec) + lam*sum(w.vec*abs(b1.vec))
                    double ob1 = obj_fun(y_vec, x_mat, b1_vec, r_eff) + arma::sum(t_vec%b1_vec) + lam*arma::sum(w_vec%arma::abs(b1_vec));
                    //ob2 = obj.fun(y.vec,x.mat,b2.vec,r.eff) + sum(t.vec*b2.vec) + lam*sum(w.vec*abs(b2.vec))
                    double ob2 = obj_fun(y_vec, x_mat, b2_vec, r_eff) + arma::sum(t_vec%b2_vec) + lam*arma::sum(w_vec%arma::abs(b2_vec));
                    //#print(c(ob1,ob2,abs(ob1-ob2)))
                    //if(ob1 < ob2){ b.vec = b2.vec } else { ob.vec = b1.vec }   #### Mod
                    if(ob1 < ob2){ b_vec = b2_vec; } else { ob_vec = b1_vec; }   // Mod
                    //if(sum(abs(b.vec-ob.vec))<b.eps) break
                    if(arma::sum(arma::abs(b_vec-ob_vec))<b_eps) { break; }
               }
               //Rcout << "CPP p_ncpen_fun | iiter:" << iiter << std::endl;
          } //mlqa
          //#################################################### if - for Add
          //#################################################### if - for Add

          //f.vec[iter] = obj.fun(y.vec,x.mat,b.vec,r.eff)+sum(w.vec*pen.fun(b.vec,lam,gam,tau))    ### (Mod) r.eff
          f_vec(iter) = obj_fun(y_vec, x_mat, b_vec, r_eff) + arma::sum(w_vec%pen_fun(b_vec, lam, gam, tau));

          //########## moved to inside for loop; was outside
          //g.vec = obj.grad.fun(y.vec,x.mat,b.vec,r.eff)+w.vec*pen.grad.fun(b.vec,lam,gam,tau)-lam*w.vec*sign(b.vec)  ### Mod: r.eff
          g_vec = obj_grad_fun(y_vec, x_mat, b_vec, r_eff) + w_vec%pen_grad_fun(b_vec, lam, gam, tau) - lam*w_vec%arma::sign(b_vec);
          //p.vec = lam*sign(w.vec*b.vec)
          //p_vec = lam*arma::sign(w_vec%b_vec);
          //p.vec = lam*w.vec*sign(b.vec) ##### sign(w.vec*b.vec) => lam*w.vec*sign(b.vec)
          p_vec = lam*w_vec%arma::sign(b_vec);

          //a.set = c(1:p)[b.vec!=0]                                            # Add
          arma::uvec a_set = arma::find(b_vec!=0);
          //n.set = c(1:p)[b.vec==0]                                            # Add
          arma::uvec n_set = arma::find(b_vec==0);
          //kkt0 = sum(abs(g.vec[a.set]+p.vec[a.set])<k.eps)==length(a.set)     # Mod sum(abs( )<k.eps)==  a.set; last sum => length
          //kkt1 = sum(abs(g.vec[n.set])-lam*w.vec[n.set]<k.eps)==length(n.set) #### (Mod) check parentheses. last sum => length
          kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps) == a_set.n_rows;
          kkt1 = arma::sum(arma::abs(g_vec(n_set))-lam*w_vec(n_set)<k_eps)== n_set.n_rows;
          //kkt0 = get_kkt0(b_vec, g_vec, p_vec, k_eps);
          //kkt1 = get_kkt1(b_vec, w_vec, g_vec, p_vec, lam, k_eps);
          //########## moved to inside for loop; was outside

          //#if(sum(abs(b.vec-ob.vec))<b.eps) break                            ### (Del)
          //if(kkt0&kkt1) break #iter break                                    ### (Add)
          if(kkt0&&kkt1) {
               iter++; // for f_vec.head(iter);
               break; // iter break                                    ### (Add)
          }
          //ob.vec = b.vec
          ob_vec = b_vec;
     }
     //Rcout << "CPP p_ncpen_fun | iter:" << iter << std::endl;

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

// ncpen.fun = function(y.vec,
//                      x.mat,x.std,intc,
//                      w.vec,
//                      lam.vec, r.lam,
//                      gam,tau,
//                      p.max,iter.max,b.eps,k.eps,p.eff,r.eff){

// // [[Rcpp::export]]
// Rcpp::List ncpen_fun(arma::vec& y_vec, arma::mat& x_mat0, bool x_std, bool intc,
//                                  arma::vec& w_vec0, arma::vec& lam_vec0, double r_lam,
//                                  double gam, double tau, arma::uword p_max, double iter_max, double b_eps, double k_eps,
//                                  arma::uword p_eff, double r_eff,
//                                  SEXP family, SEXP penalty) {
int ncpen_fun(arma::vec& y_vec, arma::mat& x_mat0, bool x_std, bool intc,
                                 arma::vec& w_vec0, arma::vec& lam_vec0, double r_lam,
                                 double gam, double tau, arma::uword p_max, double iter_max, double b_eps, double k_eps,
                                 arma::uword p_eff, double r_eff,
                                 std::string fam, std::string pen, ncpen_ret& ret_buff) {

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

     //n = dim(x.mat)[1]
     arma::uword n = x_mat.n_rows;
     //p = dim(x.mat)[2]
     arma::uword p = x_mat.n_cols;
     //r = length(lam.vec)
     arma::uword r = lam_vec.n_rows;

     // if(x.std==TRUE){
     //      std = sqrt(colSums(x.mat^2)/n)
     //      std[std==0] = 1
     //      x.mat = sweep(x.mat,2,std,"/")
     // } else {
     //     std = rep(1,p)
     // }

     arma::vec std = arma::ones<arma::vec>(p); // Should be col vec
     if(x_std == true) {
          std = arma::trans(arma::sqrt(arma::sum(arma::square(x_mat), 0)/n)); // trans row bec to col vec
          std(arma::find(std==0)).fill(1);
          x_mat = x_mat.each_row() / arma::trans(std); // trans col vec to row vec
     }
     // } else {
     //      std = rep(1,p)
     // }

     // if(intc==TRUE){
     if(intc==true) {
     //      x.mat = cbind(1,x.mat)
          x_mat = arma::join_horiz(arma::ones<arma::vec>(n), x_mat);
     //      w.vec = c(0,w.vec)
          w_vec = arma::join_vert(arma::zeros<arma::vec>(1), w_vec);
     //      std = c(1,std)
          std = arma::join_vert(arma::ones<arma::vec>(1), std);
     //      p = p+1
          p = p+1;
     }

     // if(lam.vec[1]==-1){
     if(lam_vec(0) < 0) {
     //      b.vec = rep(0,p)
          arma::vec b_vec = arma::zeros<arma::vec>(p);
     //      g.vec = abs(obj.grad.fun(y.vec,x.mat,b.vec,r.eff))
     // if(sum(w.vec==0)>0){
          if(arma::sum(w_vec==0)>0) {
               //      ax.mat = x.mat[,w.vec==0,drop=F]
               arma::mat ax_mat = x_mat.cols(arma::find(w_vec == 0));
               //      ab.vec = nr.fun(y.vec,ax.mat,r.eff,iter.max,b.eps)
               arma::vec ab_vec = nr_fun(fam, y_vec, ax_mat, r_eff, iter_max, b_eps);
               //      b.vec[w.vec==0] = ab.vec
               b_vec(arma::find(w_vec==0)) = ab_vec;
          }

          arma::vec g_vec = arma::abs(get_obj_grad_fun_ptr(fam)(y_vec, x_mat, b_vec, r_eff));
     //      lam.max = max(g.vec[w.vec!=0]/w.vec[w.vec!=0])
          arma::uvec idx = arma::find(w_vec!=0);
          double lam_max = arma::max(g_vec(idx)/w_vec(idx));
     //      lam.vec = exp(seq(log(lam.max),log(lam.max*r.lam),length.out=r))
          lam_vec = arma::exp(arma::linspace(std::log(lam_max), std::log(lam_max*r_lam), r));
     }

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

     //w.vec = (p-sum(w.vec==0))*w.vec/sum(w.vec)  ### (Add)
     w_vec = (double)(p-arma::sum(w_vec == 0))*w_vec/arma::sum(w_vec); // (Add)
     //b.vec = w.vec*0
     arma::vec b_vec = w_vec*0;

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
          //cat("lam=",lam,"\n")
          //Rcout << "lam=" << lam << std::endl;

          //[s 20170111]
          // //g.vec = obj.grad.fun(y.vec,x.mat,b.vec,r.eff)+w.vec*pen.grad.fun(b.vec,lam,gam,tau)-lam*sign(w.vec*b.vec)
          // arma::vec g_vec = obj_grad_fun(y_vec, x_mat, b_vec, r_eff) + w_vec%pen_grad_fun(b_vec, lam, gam, tau) - lam*arma::sign(w_vec%b_vec); // r.eff Add
          // //p.vec = lam*sign(w.vec*b.vec)
          // arma::vec p_vec = lam*arma::sign(w_vec%b_vec);

          //g.vec = obj.grad.fun(y.vec,x.mat,b.vec,r.eff)+w.vec*pen.grad.fun(b.vec,lam,gam,tau)-lam*w.vec*sign(b.vec)  ##### sign(w.vec*b.vec) => lam*w.vec*sign(b.vec)
          arma::vec g_vec = obj_grad_fun(y_vec, x_mat, b_vec, r_eff) + w_vec%pen_grad_fun(b_vec, lam, gam, tau) - lam*w_vec%arma::sign(b_vec);
          //p.vec = lam*w.vec*sign(b.vec)  ##### sign(w.vec*b.vec) => lam*w.vec*sign(b.vec)
          arma::vec p_vec = lam*w_vec%arma::sign(b_vec);
          //[e 20170111]

          //a.set = c(1:p)[b.vec!=0]
          arma::uvec a_set = arma::find(b_vec!=0);
          //n.set = c(1:p)[b.vec==0]
          arma::uvec n_set = arma::find(b_vec==0);
          //kkt0 = sum(abs(g.vec[a.set]+p.vec[a.set])<k.eps)==length(a.set)     # Mod sum(abs( )<k.eps)==  a.set; last sum => length
          //kkt1 = sum(abs(g.vec[n.set])-lam*w.vec[n.set]<k.eps)==length(n.set) #### (Mod) check parentheses. last sum => length
          kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps) == a_set.n_rows;
          kkt1 = arma::sum(arma::abs(g_vec(n_set))-lam*w_vec(n_set)<k_eps)== n_set.n_rows;
          //kkt0 = get_kkt0(b_vec, g_vec, p_vec, k_eps);
          //kkt1 = get_kkt1(b_vec, w_vec, g_vec, p_vec, lam, k_eps);
          //if(!(kkt0&kkt1)){
          if(!(kkt0&&kkt1)){
               //for(iter in 1:iter.max){#iter
               arma::uword iter = 0;
               for(; iter < iter_max; iter++) { //#iter
                    //p.set = n.set[(abs(g.vec[n.set])-lam*w.vec[n.set])>k.eps]
                    arma::uvec p_set = n_set( arma::find( (arma::abs(g_vec(n_set))-lam*w_vec(n_set))>k_eps ) );
                    //Rcout << "1" << std::endl;
                    //add = p.set[which.max(abs(g.vec[p.set]))]
                    arma::vec abs_g_p = arma::abs(g_vec(p_set));
                    if(abs_g_p.n_rows > 0) {
                         //Rcout << "abs_g_p.n_rows = " << abs_g_p.n_rows << std::endl;
                         arma::uword add = p_set(which_max(abs_g_p));
                         //Rcout << "3" << std::endl;
                         //a.set = c(a.set,add)
                         //a_set.insert_rows(a_set.n_rows, 1);
                         a_set.resize(a_set.n_rows + 1);
                         //Rcout << "4" << std::endl;
                         a_set(a_set.n_rows-1) = add;
                    }

                    //#############################################################################################################
                    //### R reconize 1 column matrix as vector.
                    // ax.mat = matrix(x.mat[,a.set],ncol=length(a.set))
                    // b.vec[a.set]  = p.ncpen.fun(y.vec,ax.mat,b.vec[a.set],w.vec[a.set],lam,gam,tau,iter.max,b.eps,k.eps,p.eff,r.eff)$b.vec
                    // cpp
                    // b.vec[a.set] = p.ncpen.fun(y.vec,x.mat[,a.set],b.vec[a.set],w.vec[a.set],lam,gam,tau,iter.max,b.eps,k.eps,p.eff,r.eff)$b.vec
                    arma::mat ax_mat = x_mat.cols(a_set);
                    arma::vec ab_vec = b_vec(a_set);
                    arma::vec aw_vec = w_vec(a_set);
                    // b_vec(a_set) = Rcpp::as<arma::vec>(
                    //      p_ncpen_fun(y_vec, ax_mat, ab_vec, aw_vec, lam, gam, tau, iter_max, b_eps, k_eps, p_eff, r_eff,
                    //                  obj_fun, obj_grad_fun, obj_hess_fun, pen_fun, pen_grad_fun)["b.vec"]);
                    //b_vec(a_set).print();

                    p_ncpen_ret ret_buff;
                    p_ncpen_fun(y_vec, ax_mat, ab_vec, aw_vec, lam, gam, tau, iter_max, b_eps, k_eps, p_eff, r_eff,
                                     obj_fun, obj_grad_fun, obj_hess_fun, pen_fun, pen_grad_fun, ret_buff);
                    b_vec(a_set) = ret_buff.b_vec;

                    //#############################################################################################################

                    //g.vec = obj.grad.fun(y.vec,x.mat,b.vec,r.eff)+w.vec*pen.grad.fun(b.vec,lam,gam,tau)-lam*sign(w.vec*b.vec)  ###r.eff Add
                    g_vec = obj_grad_fun(y_vec, x_mat, b_vec, r_eff) + w_vec%pen_grad_fun(b_vec, lam, gam, tau) - lam*w_vec%arma::sign(b_vec);
                    //p.vec = lam*sign(w.vec*b.vec)
                    p_vec = lam*w_vec%arma::sign(b_vec);
                    //a.set = c(1:p)[b.vec!=0]
                    a_set = arma::find(b_vec!=0);
                    //n.set = c(1:p)[b.vec==0]
                    n_set = arma::find(b_vec==0);
                    //kkt0 = sum(abs(g.vec[a.set]+p.vec[a.set])<k.eps)==length(a.set)     # Mod sum(abs( )<k.eps)==  a.set; last sum => length
                    //kkt1 = sum(abs(g.vec[n.set])-lam*w.vec[n.set]<k.eps)==length(n.set) #### (Mod) check parentheses. last sum => length
                    kkt0 = arma::sum(arma::abs(g_vec(a_set)+p_vec(a_set))<k_eps) == a_set.n_rows;
                    kkt1 = arma::sum(arma::abs(g_vec(n_set))-lam*w_vec(n_set)<k_eps)== n_set.n_rows;
                    //kkt0 = get_kkt0(b_vec, g_vec, p_vec, k_eps);
                    //kkt1 = get_kkt1(b_vec, w_vec, g_vec, p_vec, lam, k_eps);
                    //if(kkt0&kkt1) break #iter break
                    if(kkt0&&kkt1) { break; } // iter break
               } //iter
               // Rcout <<"ncpen_fun | pos:" <<pos+1<<"    iter:"<<iter+1<<endl;
          }

          //b.mat[,pos] = b.vec
          b_mat.col(pos) = b_vec;
          //g.mat[,pos] = g.vec+p.vec
          g_mat.col(pos) = g_vec + p_vec;
          //c.mat[,pos] = c(kkt0,kkt1)
          c_mat(0, pos) = kkt0; c_mat(1, pos) = kkt1;
          //f.vec[pos] = obj.fun(y.vec,x.mat,b.vec,r.eff)+sum(w.vec*pen.fun(b.vec,lam,gam,tau))     ###r.eff Add
          f_vec(pos) = obj_fun(y_vec, x_mat, b_vec, r_eff)+arma::sum(w_vec%pen_fun(b_vec, lam, gam, tau));
          //d.vec[pos] = length(a.set)    ### (Mod) sum => length
          a_set = arma::find(b_vec!=0);
          d_vec(pos) = a_set.n_rows;

          // Mod
          //if(length(a.set)>p.max){
          if(a_set.n_rows > p_max) {
               //cat("(warning) ncpen.fun stops because the number of parameters exceeds",p.max,"\n")
               //Rcout << "(Warning) ncpen_fun stops because the number of parameters exceeds " << p_max << "." << std::endl;
               //cat("(warning) increase p.max", "\n")                                        ### (Mod)
               //Rcout << "(Warning) increase p_max." << std::endl;
               warning_code_vec(0) = 1;
               //break #pos break
               pos++; // for head_cols functions
               break; // pos break
          }

          // Add
          // if(fam=="log"){
          if(obj_fun == log_obj_fun) {
               //if(max(abs(b.vec))>5){
               if(arma::max(arma::abs(b_vec))>50){
                    //cat("(warning) ncpen.fun stops because a paramter diverges", "\n")
                    //Rcout << "(Warning) ncpen_fun stops because a paramter diverges." << std::endl;
                    //cat("(warning) increase r.eff", "\n")
                    //Rcout << "(Warning) increase r_eff." << std::endl;
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
     //Rcout <<"pos:" << pos << "   b_mat.n_cols:" << b_mat.n_cols << "   b_mat.n_rows:" <<b_mat.n_rows <<endl;
     if(x_std==true){
          //        Rcout << b_mat <<endl;
          b_mat = b_mat.each_col()/std;
          //Rcout << b_mat <<endl;
          //Rcout << std <<endl;

     }

     arma::mat ret_b_mat = b_mat.head_cols(pos);
     arma::mat ret_g_mat = g_mat.head_cols(pos);
     arma::vec ret_f_vec = f_vec.head(pos);
     arma::mat ret_c_mat = c_mat.head_cols(pos);
     arma::vec ret_lam_vec = lam_vec.head(pos);
     arma::vec ret_d_vec = d_vec.head(pos);

     // return Rcpp::List::create(Rcpp::Named("b.mat") = ret_b_mat,
     //                           Rcpp::Named("g.mat") = ret_g_mat,
     //                           Rcpp::Named("f.vec") = ret_f_vec,
     //                           Rcpp::Named("c.mat") = ret_c_mat,
     //                           Rcpp::Named("lam.vec") = ret_lam_vec,
     //                           Rcpp::Named("d.vec") = ret_d_vec,
     //                           Rcpp::Named("w.vec") = w_vec,
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
