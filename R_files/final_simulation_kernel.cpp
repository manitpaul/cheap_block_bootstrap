#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Store mass, first moments and upper-triangular second moments.
// [[Rcpp::export]]
NumericMatrix precision_features(NumericMatrix x) {
 int n=x.nrow(),p=x.ncol(),q=1+p+p*(p+1)/2;NumericMatrix f(n,q);
 for(int i=0;i<n;i++){f(i,0)=1;for(int j=0;j<p;j++)f(i,1+j)=x(i,j);
 int k=1+p;for(int j=0;j<p;j++)for(int h=0;h<=j;h++)f(i,k++)=x(i,j)*x(i,h);}
 return f;
}

// Exact full/weighted statistic for each population scale alternative.
std::vector<double> transform(const std::vector<double>& u,int kind,int p,
                             const NumericVector& alternatives) {
 int d=alternatives.size();std::vector<double> out(d);
 if(kind<2){double mean=u[1]/u[0];for(int k=0;k<d;k++)out[k]=kind==0?std::log(mean)+alternatives[k]:std::min(1000.0,alternatives[k]/mean);return out;}
 arma::mat cov(p,p);std::vector<double> mu(p);for(int j=0;j<p;j++)mu[j]=u[1+j]/u[0];
 int z=1+p;for(int j=0;j<p;j++)for(int h=0;h<=j;h++)cov(j,h)=cov(h,j)=u[z++]/u[0]-mu[j]*mu[h];
 // Only the FIRST population coordinate changes scale. The Schur complement
 // makes all alternatives exact after one nuisance-block Cholesky factorization.
 int m=p-1;std::vector<double>L(m*m),v(m);
 for(int i=0;i<m;i++)for(int j=0;j<=i;j++){double z=cov(i+1,j+1)+(i==j?.1:0);
 for(int h=0;h<j;h++)z-=L[i*m+h]*L[j*m+h];
 if(i==j){if(!(z>0))stop("Nonpositive regularized covariance");L[i*m+j]=std::sqrt(z);}else L[i*m+j]=z/L[j*m+j];}
 double q=cov(0,0);for(int i=0;i<m;i++){double z=cov(i+1,0);for(int j=0;j<i;j++)z-=L[i*m+j]*v[j];v[i]=z/L[i*m+i];q-=v[i]*v[i];}
 if(q < -1e-9)stop("Negative covariance Schur complement");q=std::max(q,0.0);
 for(int k=0;k<d;k++)out[k]=1/(alternatives[k]*q+.1);
 return out;
}

// Draw indices for ordinary/m-out-of-n/CBB summaries. All mass is normalized.
// [[Rcpp::export]]
NumericMatrix draw_features(NumericMatrix f,int draws,int B,int kind,int p,NumericVector alternatives) {
 RNGScope rng;int n=f.nrow(),q=f.ncol(),d=alternatives.size();NumericMatrix ans(B,d);
 std::vector<double>a((size_t)n*q),u(q);for(int i=0;i<n;i++)for(int j=0;j<q;j++)a[(size_t)i*q+j]=f(i,j);
 for(int b=0;b<B;b++){std::fill(u.begin(),u.end(),0.0);for(int k=0;k<draws;k++){int i=R::unif_rand()*n;for(int j=0;j<q;j++)u[j]+=a[(size_t)i*q+j];}
 auto z=transform(u,kind,p,alternatives);for(int j=0;j<d;j++)ans(b,j)=z[j];if(b%100==0)checkUserInterrupt();}return ans;
}

// [[Rcpp::export]]
NumericVector feature_estimate(NumericMatrix f,int kind,int p,NumericVector alternatives) {
 std::vector<double>u(f.ncol());for(int i=0;i<f.nrow();i++)for(int j=0;j<f.ncol();j++)u[j]+=f(i,j);
 return wrap(transform(u,kind,p,alternatives));
}

// BLB counts have TOTAL n, not b. Consecutive input rows form disjoint bags.
// Results are already centered at each bag's plugin parameter.
// [[Rcpp::export]]
List draw_blb_features(NumericMatrix f,IntegerVector sizes,int n,int B,int kind,int p,NumericVector alternatives) {
 RNGScope rng;int q=f.ncol(),d=alternatives.size(),start=0;List ans(sizes.size());
 for(int k=0;k<sizes.size();k++){int bsize=sizes[k];std::vector<double>a((size_t)bsize*q),u(q),v(q);
 for(int i=0;i<bsize;i++)for(int j=0;j<q;j++){a[(size_t)i*q+j]=f(start+i,j);v[j]+=f(start+i,j);}
 auto center=transform(v,kind,p,alternatives);NumericMatrix out(B,d);
 for(int r=0;r<B;r++){std::fill(u.begin(),u.end(),0.0);int left=n;
 for(int i=0;i<bsize;i++){int w=i==bsize-1?left:(int)R::rbinom(left,1.0/(bsize-i));left-=w;for(int j=0;j<q;j++)u[j]+=w*a[(size_t)i*q+j];}
 auto z=transform(u,kind,p,alternatives);for(int j=0;j<d;j++)out(r,j)=z[j]-center[j];}
 ans[k]=out;start+=bsize;if(k%50==0)checkUserInterrupt();}return ans;
}

// Same BLB experiment as above, accumulating the desired quality assessments
// directly to avoid materializing millions of redundant resample values in R.
// [[Rcpp::export]]
NumericMatrix blb_feature_summary(NumericMatrix f,IntegerVector sizes,int n,int B,int kind,int p,NumericVector alternatives) {
 RNGScope rng;int q=f.ncol(),d=alternatives.size(),start=0,K=sizes.size();NumericMatrix ans(3,d);
 int low=std::floor((B+1)*.025+1e-10)-1,high=low+std::ceil((B+1)*.95-1e-10);
 for(int k=0;k<K;k++){int bs=sizes[k];std::vector<double>a((size_t)bs*q),u(q),v(q),values((size_t)B*d);
 for(int i=0;i<bs;i++)for(int j=0;j<q;j++){a[(size_t)i*q+j]=f(start+i,j);v[j]+=f(start+i,j);}auto center=transform(v,kind,p,alternatives);
 for(int r=0;r<B;r++){std::fill(u.begin(),u.end(),0.0);int left=n;for(int i=0;i<bs;i++){int w=i==bs-1?left:(int)R::rbinom(left,1.0/(bs-i));left-=w;for(int j=0;j<q;j++)u[j]+=w*a[(size_t)i*q+j];}
 auto z=transform(u,kind,p,alternatives);for(int j=0;j<d;j++)values[(size_t)j*B+r]=z[j]-center[j];}
 for(int j=0;j<d;j++){auto begin=values.begin()+(size_t)j*B;double mean=0,m2=0;for(int r=0;r<B;r++){double e=begin[r]-mean;mean+=e/(r+1);m2+=e*(begin[r]-mean);}ans(0,j)+=std::sqrt(m2/(B-1))/K;
 std::nth_element(begin,begin+low,begin+B);ans(1,j)+=begin[low]/K;std::nth_element(begin,begin+high,begin+B);ans(2,j)+=begin[high]/K;}
 start+=bs;if(k%100==0)checkUserInterrupt();}return ans;
}

// Exact Euclidean-distance U kernel. No finite-moment reduction is used.
// [[Rcpp::export]]
NumericMatrix distance_kernel(NumericMatrix x) {
 int n=x.nrow(),p=x.ncol();NumericMatrix h(n,n);for(int i=0;i<n;i++)for(int j=0;j<i;j++){
 double z=0;for(int k=0;k<p;k++)z+=std::pow(x(i,k)-x(j,k),2);h(i,j)=h(j,i)=std::sqrt(z);}return h;
}

// [[Rcpp::export]]
NumericMatrix compress_kernel(NumericMatrix h,IntegerVector groups,int s) {
 int n=h.nrow();NumericMatrix a(s,s);for(int i=0;i<n;i++)for(int j=0;j<n;j++)a(groups[i],groups[j])+=h(i,j);return a;
}

// Streaming preprocessing demonstrates the advertised CBB storage O(s^2):
// no n by n matrix is retained. With s=0, return only the common point fit.
// [[Rcpp::export]]
List distance_bin_cache(NumericMatrix x,IntegerVector groups,int s) {
 int n=x.nrow(),p=x.ncol();NumericMatrix a(s,s);double total=0;
 for(int i=0;i<n;i++)for(int j=0;j<i;j++){double z=0;for(int k=0;k<p;k++)z+=std::pow(x(i,k)-x(j,k),2);z=std::sqrt(z);total+=2*z;
 if(s>0){a(groups[i],groups[j])+=z;a(groups[j],groups[i])+=z;}}
 return List::create(_["table"]=a,_["estimate"]=total/(n*(n-1.0)));
}

// [[Rcpp::export]]
NumericVector draw_u(NumericMatrix h,IntegerVector masses,int draws,int B) {
 RNGScope rng;int s=h.nrow();NumericVector out(B);arma::mat H(h.begin(),s,s,false);
 for(int first=0;first<B;first+=100){int nr=std::min(100,B-first);arma::mat W(s,nr,arma::fill::zeros);
 for(int r=0;r<nr;r++)for(int t=0;t<draws;t++)W((int)(R::unif_rand()*s),r)++;
 arma::mat HW=H*W;
 for(int r=0;r<nr;r++){double total=0,sum=0;for(int i=0;i<s;i++){total+=W(i,r)*masses[i];sum+=W(i,r)*HW(i,r);}out[first+r]=sum/(total*(total-1));}
 checkUserInterrupt();}return out;
}

// m/n U draws compute only the selected m by m kernel, preserving fairness.
// [[Rcpp::export]]
NumericVector draw_u_small(NumericMatrix x,int draws,int B) {
 RNGScope rng;int n=x.nrow(),p=x.ncol();NumericVector out(B);std::vector<int>ix(draws);
 for(int r=0;r<B;r++){for(int t=0;t<draws;t++)ix[t]=R::unif_rand()*n;double sum=0;
 for(int i=0;i<draws;i++)for(int j=0;j<i;j++){double z=0;for(int k=0;k<p;k++)z+=std::pow(x(ix[i],k)-x(ix[j],k),2);sum+=std::sqrt(z);}out[r]=2*sum/(draws*(draws-1.0));}return out;
}

// [[Rcpp::export]]
List draw_blb_u(NumericMatrix x,IntegerVector sizes,int n,int B) {
 RNGScope rng;int start=0,p=x.ncol();List ans(sizes.size());
 for(int k=0;k<sizes.size();k++){int b=sizes[k];std::vector<double>h(b*b);double plugin=0;
 for(int i=0;i<b;i++)for(int j=0;j<i;j++){double z=0;for(int v=0;v<p;v++)z+=std::pow(x(start+i,v)-x(start+j,v),2);h[i*b+j]=h[j*b+i]=std::sqrt(z);plugin+=2*std::sqrt(z);}plugin/=b*(double)b;
 NumericVector out(B);std::vector<int>w(b);
 for(int r=0;r<B;r++){int left=n;for(int i=0;i<b;i++){w[i]=i==b-1?left:(int)R::rbinom(left,1.0/(b-i));left-=w[i];}
 double sum=0;for(int i=0;i<b;i++)for(int j=0;j<i;j++)sum+=2.0*w[i]*w[j]*h[i*b+j];out[r]=sum/(n*(n-1.0))-plugin;}
 ans[k]=out;start+=b;}return ans;
}

